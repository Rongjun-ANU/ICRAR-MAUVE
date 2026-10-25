"""Approximate mosaic geometry and regional spectral diagnostics (exploratory)."""
from pathlib import Path
import json, warnings, re
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import ConvexHull
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks
from astropy.io import fits
from astropy.wcs import WCS
import hashlib, csv

AO_GAP=(5770.,6050.)
WINDOWS=(5100.,5500.,6400.,6800.)

def allowed_wavelength(wave, ao=False):
    wave=np.asarray(wave)
    return np.isfinite(wave) & (~((wave>=AO_GAP[0])&(wave<=AO_GAP[1])) if ao else True)

def wavelength(header,n):
    if header.get('CTYPE3','').strip() not in ('AWAV','WAVE','LINEAR',''):
        raise ValueError('Unsupported spectral WCS')
    unit=header.get('CUNIT3','Angstrom').strip().lower()
    factor={'angstrom':1,'ang':1,'nm':10,'m':1e10}.get(unit)
    if factor is None: raise ValueError('Unsupported wavelength unit: '+unit)
    return (header['CRVAL3']+(np.arange(n)+1-header['CRPIX3'])*header.get('CD3_3',header.get('CDELT3')))*factor

def normalized_ratio(wave,a,b,norm_band=(5400,5600),ao=False):
    wave=np.asarray(wave); a=np.asarray(a); b=np.asarray(b)
    keep=allowed_wavelength(wave,ao)&np.isfinite(a)&np.isfinite(b)&(a>0)&(b>0)
    norm=keep&(wave>=norm_band[0])&(wave<=norm_band[1])
    r=np.full(a.shape,np.nan)
    if not np.any(norm): return r
    an=np.median(a[norm]); bn=np.median(b[norm])
    if an>0 and bn>0: r[keep]=(a[keep]/an)/(b[keep]/bn)
    return r

def rotate(x,y,angle):
    t=np.deg2rad(angle)
    return x*np.cos(t)+y*np.sin(t),-x*np.sin(t)+y*np.cos(t)

def _profile(z,coord,mask):
    idx=np.floor(coord[mask]).astype(int)
    lo,hi=idx.min(),idx.max()
    order=np.argsort(idx); idx=idx[order]; vals=z[mask][order]
    cuts=np.r_[0,np.where(np.diff(idx))[0]+1,len(idx)]
    p=np.full(hi-lo+1,np.nan)
    for i,j in zip(cuts[:-1],cuts[1:]):
        if j-i>=10: p[idx[i]-lo]=np.median(vals[i:j])
    good=np.isfinite(p)
    p=np.interp(np.arange(len(p)),np.where(good)[0],p[good])
    return np.arange(lo,hi+1)+.5,ndi.gaussian_filter1d(p,2)

def segment_variance(maps,pixel_scale=.2,prominence=.055):
    """Rectangular common-orientation candidate grid; never claimed as OB truth.

    Missing cells are allowed. Boundary candidates must recur in >=3 windows,
    and have tile spacing compatible with approximately 60 arcsec footprints.
    """
    stack=np.array(maps,dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        z=np.nanmedian(np.log10(np.where(stack>0,stack,np.nan)),axis=0)
    valid=np.isfinite(z)
    if valid.sum()<100: raise ValueError('Insufficient valid variance pixels')
    y,x=np.indices(z.shape)
    pts=np.column_stack([x[valid],y[valid]])
    hull=pts[ConvexHull(pts).vertices]
    def area(a):
        u,v=rotate(hull[:,0],hull[:,1],a)
        return np.ptp(u)*np.ptp(v)
    angles=np.arange(-45,45,.5)
    a0=angles[np.argmin([area(a) for a in angles])]
    angle=minimize_scalar(area,bounds=(a0-1,a0+1),method='bounded').x
    u,v=rotate(x,y,angle)
    bounds=[(float(q[valid].min()),float(q[valid].max())) for q in [u,v]]
    side=60/pixel_scale
    cuts=[]; diagnostics=[]
    for axis,q in enumerate([u,v]):
        profs=[]
        # Exclude external rim and suppress narrow single-spaxel defects.
        inner=ndi.binary_erosion(valid,iterations=max(2,int(2/pixel_scale)))
        for m in stack:
            zz=np.log10(np.where(m>0,m,np.nan))
            coords,p=_profile(zz,q,inner&np.isfinite(zz))
            trend=ndi.gaussian_filter1d(p,max(10,side/8))
            res=p-trend
            # Both lower and higher variance seams occur in real mosaics.
            strength=np.abs(res)
            peaks,props=find_peaks(strength,prominence=prominence,distance=max(15,int(side*.35)))
            peaks=[k for k in peaks if coords[k]>bounds[axis][0]+side*.55 and coords[k]<bounds[axis][1]-side*.55]
            profs.append((coords,p,res,[float(coords[k]) for k in peaks]))
        candidates=sorted({c for pr in profs for c in pr[3]})
        groups=[]
        tolerance=side*.12
        for c in candidates:
            if not groups or c-np.median(groups[-1])>tolerance: groups.append([c])
            else: groups[-1].append(c)
        accepted=[]
        for g in groups:
            center=float(np.median(g))
            votes=sum(any(abs(c-center)<=tolerance for c in pr[3]) for pr in profs)
            if votes>=min(3,len(profs)): accepted.append(center)
        # Reject closely spaced spurious features; retain stronger variance evidence.
        accepted2=[]
        for c in accepted:
            if not accepted2 or c-accepted2[-1]>=side*.6: accepted2.append(c)
        cuts.append(accepted2)
        diagnostics.append(profs)
    labels=np.zeros(z.shape,dtype='i2'); regions=[]; cells=[]
    margin=max(8,int(4/pixel_scale))
    edges=[[b[0],*c,b[1]] for b,c in zip(bounds,cuts)]
    local_v_cuts={}
    rectangles=[]
    for i,(u0,u1) in enumerate(zip(edges[0][:-1],edges[0][1:])):
        # A staggered mosaic can have different horizontal joins per column.
        select=valid&(u>u0+margin)&(u<u1-margin)
        local=[]
        if select.sum()>100:
            xx,pp=_profile(z,v,select)
            rr=pp-ndi.gaussian_filter1d(pp,max(10,side/8))
            peaks,_=find_peaks(abs(rr),prominence=prominence,distance=max(15,int(side*.35)))
            possible=[float(xx[k]) for k in peaks if xx[k]>xx.min()+side*.5 and xx[k]<xx.max()-side*.5]
            for cut in cuts[1]:
                near=[c for c in possible if abs(c-cut)<side*.42]
                if near: local.append(min(near,key=lambda c:abs(c-cut)))
            local=sorted(set(local))
        local_v_cuts[i]=local
    # Mask empty cells of irregular mosaics, rather than force a full rectangle.
    for i,(u0,u1) in enumerate(zip(edges[0][:-1],edges[0][1:])):
        ve=[bounds[1][0],*local_v_cuts[i],bounds[1][1]]
        for j,(v0,v1) in enumerate(zip(ve[:-1],ve[1:])):
            cell=(u>=u0)&(u<u1)&(v>=v0)&(v<v1)
            interior=valid&(u>u0+margin)&(u<u1-margin)&(v>v0+margin)&(v<v1-margin)
            occupancy=(valid&cell).sum()/max(1,cell.sum())
            if interior.sum()<side**2*.2 or occupancy<.4: continue
            regions.append(interior); cells.append((i,j));rectangles.append((u0,u1,v0,v1));labels[interior]=len(regions)
    # Shared footprints expand each geometric cell into its neighbours. The
    # 4 arcsec half-width is explicit and adjustable, not exact exposure support.
    footprints=[]
    for u0,u1,v0,v1 in rectangles:
        footprints.append(valid&(u>=u0-margin)&(u<=u1+margin)&(v>=v0-margin)&(v<=v1+margin))
    adjacent=[]
    for a in range(len(cells)):
        for b in range(a+1,len(cells)):
            aa=rectangles[a];bb=rectangles[b]
            du=min(aa[1],bb[1])-max(aa[0],bb[0]);dv=min(aa[3],bb[3])-max(aa[2],bb[2])
            # Require a shared edge of meaningful length; reject corner-only joins.
            if max(du,dv)>side*.25 and min(du,dv)>-1 and np.any(footprints[a]&footprints[b]): adjacent.append((a,b))
    widths=[np.diff(e) for e in edges]
    size_ok=all(np.all((ww>side*.5)&(ww<side*1.5)) for ww in widths)
    n=len(regions)
    status='single_candidate' if n==1 and all(b[1]-b[0]<side*1.5 for b in bounds) else ('multiple_candidate' if n>1 and size_ok else 'review_required')
    return dict(n_regions=n,status=status,angle=float(angle),bounds=bounds,cuts=cuts,cells=cells,
                labels=labels,footprints=np.asarray(footprints,dtype='u1'),interiors=np.asarray(regions),
                adjacent=adjacent,variance_log=z,profiles=diagnostics,pixel_scale=pixel_scale,
                overlap_halfwidth_pixels=margin,local_v_cuts=local_v_cuts,rectangles=rectangles)

def sample_cube(path,cache=None):
    path=Path(path)
    if cache and Path(cache).exists():
        with np.load(cache) as d: return {k:d[k] for k in d.files}
    out={}
    with fits.open(path,memmap=True) as h:
        hd=h['DATA'].header;wave=wavelength(hd,h['DATA'].shape[0])
        out['header']=np.array(hd.tostring());out['primary']=np.array(h[0].header.tostring())
        out['ao']=np.array('AO' in str(h[0].header.get('ESO INS MODE','')))
        for c in WINDOWS:
            idx=np.where((abs(wave-c)<=5)&allowed_wavelength(wave,bool(out['ao'])))[0]
            if len(idx)<3: raise ValueError(f'Window {c} unavailable')
            f=np.array(h['DATA'].data[idx],dtype=float);v=np.array(h['STAT'].data[idx],dtype=float)
            good=np.isfinite(f)&np.isfinite(v)&(v>0)
            if 'DQ' in h: good &= h['DQ'].data[idx]==0
            v[~good]=np.nan;f[~good]=np.nan
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                out[f'var_{int(c)}']=np.nanmedian(v,axis=0).astype('f4')
                out[f'flux_{int(c)}']=np.nanmedian(f,axis=0).astype('f4')
    if cache: np.savez_compressed(cache,**out)
    return out

def discover_products(root):
    """CONFIG is authoritative for inputs and actual wavelength ranges."""
    import yaml
    records=[]
    for run,folder in [('7000','v3tk_v7.6.8_7000'),('8900','v3tk_v7.6.8')]:
        for cp in sorted((Path(root)/folder).glob('*/CONFIG')):
            try:
                cfg=yaml.safe_load(cp.read_text());gal=str(cfg['GENERAL']['RUN_ID'])
                if re.search(r'_(SFH\d+|SNR\d+)',gal): continue
                cube=Path(cfg['GENERAL']['INPUT'])
                cont=next(iter(sorted(cp.parent.glob(gal+'_cont_cube.fits*'))),None)
                records.append(dict(galaxy=gal,run=run,cube=str(cube),continuum=str(cont) if cont else None,
                                    config=str(cp),ao='AO' in str(cfg['READ_DATA']['METHOD']),
                                    lmin=cfg['CONT']['LMIN'],lmax=cfg['CONT']['LMAX'],
                                    mdeg=cfg['CONT'].get('MDEG'),sfh_mdeg=cfg['SFH'].get('MDEG')))
            except Exception as exc:
                records.append(dict(galaxy=cp.parent.name,run=run,error=str(exc),config=str(cp)))
    return records

def fingerprint(path,extra=''):
    p=Path(path);st=p.stat()
    return hashlib.sha256((str(p)+str(st.st_size)+str(st.st_mtime_ns)+extra).encode()).hexdigest()[:16]

def map_interiors(result,reference_header,target_header,target_shape):
    """Map geometry by celestial WCS; reject non-integer spatial regridding."""
    src=WCS(reference_header).celestial;dst=WCS(target_header).celestial
    yy,xx=np.indices(target_shape)
    sx,sy=src.world_to_pixel_values(*dst.pixel_to_world_values(xx,yy))
    finite=np.isfinite(sx)&np.isfinite(sy)
    if np.nanmax(abs(sx[finite]-np.rint(sx[finite])))>.1 or np.nanmax(abs(sy[finite]-np.rint(sy[finite])))>.1:
        raise ValueError('Non-integer WCS reprojection requires explicit scientific review')
    ix=np.rint(np.where(finite,sx,-1)).astype(int);iy=np.rint(np.where(finite,sy,-1)).astype(int)
    ny,nx=result['labels'].shape
    good=finite&(ix>=0)&(iy>=0)&(ix<nx)&(iy<ny)
    masks=[]
    for m in result['interiors']:
        mm=np.zeros(target_shape,bool);mm[good]=m[iy[good],ix[good]];masks.append(mm)
    return masks

def extract_integrated(path,result,reference_header,ao=False,stride=4):
    """Integrate every valid interior spaxel, with constant spatial support.

    Process one region at a time; output spectra sampled every stride channels.
    A pixel is kept only if finite at EVERY retained non-AO wavelength. No
    channel-dependent nansum area changes, interpolation or random subsampling.
    """
    if str(path).endswith('.gz'):
        raise ValueError('Compressed continuum needs an uncompressed scratch copy before bulk extraction')
    with fits.open(path,memmap=True) as h:
        hd=next(q for q in h if q.data is not None and q.data.ndim==3)
        wave=wavelength(hd.header,hd.shape[0])[::stride]
        allowed=allowed_wavelength(wave,ao)
        masks=map_interiors(result,reference_header,hd.header,hd.shape[1:])
        # Interpolated CONT cubes can have unsupported endpoint channels. A
        # deterministic spatial preflight identifies common usable channels;
        # otherwise one globally-NaN plane would reject every good spaxel.
        coverage=[]
        for m in masks:
            yy,xx=np.where(m)
            take=np.linspace(0,len(xx)-1,min(len(xx),1024),dtype=int)
            if not len(take):
                coverage.append(np.zeros(len(wave)));continue
            probe=np.asarray(hd.data[::stride,yy[take],xx[take]],dtype='f4')
            coverage.append(np.isfinite(probe).mean(axis=1))
        coverage=np.asarray(coverage)
        allowed &= np.min(coverage,axis=0)>=.5
        if allowed.sum()<3: raise ValueError('Too few commonly supported spectral channels')
        spectra=[];counts=[];totals=[]
        for m in masks:
            yy,xx=np.where(m);totals.append(len(xx))
            if len(xx)==0:
                spectra.append(np.full(len(wave),np.nan));counts.append(0);continue
            values=np.asarray(hd.data[::stride,yy,xx],dtype='f4')
            support=np.all(np.isfinite(values[allowed]),axis=0)
            n=int(support.sum());counts.append(n)
            s=np.full(len(wave),np.nan)
            if n>=100:
                s[allowed]=np.sum(values[allowed][:,support],axis=1,dtype='f8')
            spectra.append(s)
            del values
    return dict(wave=wave,spectra=np.asarray(spectra),counts=np.asarray(counts),candidate_counts=np.asarray(totals),
                supported_channels=allowed,preflight_finite_fraction=coverage)

def plot_geometry(gal,samples,result,out):
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,3,figsize=(16,10),constrained_layout=True)
    for ax,w in zip(axs.flat[:4],WINDOWS):
        z=np.log10(samples[f'var_{int(w)}']);lo,hi=np.nanpercentile(z,[5,95])
        im=ax.imshow(z,origin='lower',vmin=lo,vmax=hi,cmap='viridis')
        ax.set_title(f'{int(w)} A: log variance')
        fig.colorbar(im,ax=ax,shrink=.6)
    axs[1,1].imshow(np.ma.masked_equal(result['labels'],0),origin='lower',cmap='tab20')
    for i,m in enumerate(result['interiors']):
        yy,xx=np.where(m);axs[1,1].text(np.median(xx),np.median(yy),f'R{i+1}',ha='center')
    axs[1,1].set_title('Interiors: overlaps and rim excluded')
    members=result['footprints'].sum(axis=0)
    im=axs[1,2].imshow(np.ma.masked_equal(members,0),origin='lower',vmin=1,vmax=max(2,members.max()))
    fig.colorbar(im,ax=axs[1,2],shrink=.6,label='Shared footprint memberships')
    axs[1,2].set_title('Approximate shared footprint support')
    fig.suptitle(f'{gal}: {result["status"]}, {result["n_regions"]} candidate regions; angle={result["angle"]:.2f} deg')
    fig.savefig(out,dpi=140);plt.close(fig)

def save_masks(path,result,header):
    wh=WCS(header).celestial.to_header();hdus=[fits.PrimaryHDU()]
    wh['COMMENT']='Approximate variance geometry; R labels are not verified OB IDs.'
    hdus.append(fits.ImageHDU(result['labels'],header=wh,name='INTERIORS'))
    for i,m in enumerate(result['footprints']):
        hdus.append(fits.ImageHDU(m,header=wh,name=f'FOOTPRINT_R{i+1}'))
    hdus.append(fits.ImageHDU(result['footprints'].sum(axis=0).astype('u1'),header=wh,name='MEMBERS'))
    fits.HDUList(hdus).writeto(path,overwrite=True)

def plot_ratios(gal,result,run_spectra,records,outdir):
    import matplotlib.pyplot as plt
    pairs=result['adjacent'];rows=[]
    if not pairs: return rows
    fig,axs=plt.subplots(len(pairs),2,figsize=(14,max(4,2.25*len(pairs))),squeeze=False,constrained_layout=True)
    for col,run in enumerate(['7000','8900']):
        data=run_spectra.get(run);record=next((r for r in records if r['run']==run),{})
        for k,(a,b) in enumerate(pairs):
            ax=axs[k,col];ax.axhline(1,color='k',lw=.6)
            if data is None:
                ax.text(.5,.5,'Unavailable / failed: see status table',transform=ax.transAxes,ha='center')
                continue
            wave=data['wave'];aa=data['spectra'][a];bb=data['spectra'][b]
            ratio=normalized_ratio(wave,aa,bb,ao=record.get('ao',False))
            ax.plot(wave,ratio,lw=.7)
            if record.get('ao'): ax.axvspan(*AO_GAP,color='.8',alpha=.5)
            ax.set_title(f'{run}: R{a+1}/R{b+1}; Npix={data["counts"][a]}/{data["counts"][b]}',fontsize=10)
            ax.set_ylabel('Normalized ratio')
            finite=np.isfinite(ratio)
            if finite.any():
                q=np.percentile(ratio[finite],[1,99]);pad=max(.05,(q[1]-q[0])*.15)
                ax.set_ylim(min(q[0]-pad,.95),max(q[1]+pad,1.05))
                raw=aa/bb
                rows.append(dict(galaxy=gal,run=run,pair=f'R{a+1}/R{b+1}',n_a=int(data['counts'][a]),n_b=int(data['counts'][b]),
                                 median_raw_ratio=float(np.nanmedian(raw[finite])),ratio_p16=float(np.percentile(ratio[finite],16)),ratio_p84=float(np.percentile(ratio[finite],84))))
                np.savez_compressed(outdir/f'{gal}_{run}_R{a+1}_R{b+1}_ratio.npz',wave=wave,ratio=ratio,integrated_a=aa,integrated_b=bb)
            axs[-1,col].set_xlabel('Observed wavelength (Angstrom)')
    fig.suptitle(f'{gal}: adjacent interior integrated-continuum ratios; each spectrum normalized at 5400-5600 A\nRegional population/dust gradients remain; curves are not isolated sky-subtraction errors.',fontsize=12)
    fig.savefig(outdir/f'{gal}_ratios_7000_8900.png',dpi=130);plt.close(fig)
    return rows

def write_csv(path,rows):
    if not rows:
        Path(path).write_text('');return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)

def run_population(product_root,outdir,galaxies=None,stride=4,prominence=.055):
    """Sequential, checkpointed, read-only science inputs; errors remain visible."""
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    cache=outdir/'cache';cache.mkdir(exist_ok=True)
    records=discover_products(product_root);write_csv(outdir/'inventory.csv',records)
    statuses=[];metrics=[]
    names=sorted({r['galaxy'] for r in records}) if galaxies is None else galaxies
    for gal in names:
        rr=[r for r in records if r['galaxy']==gal and 'error' not in r]
        print('START',gal,flush=True)
        try:
            if not rr: raise ValueError('No readable CONFIG')
            rec=next((r for r in rr if r['run']=='7000'),rr[0]);cube=rec['cube']
            key=fingerprint(cube,str(WINDOWS)+'v1')
            sample=sample_cube(cube,cache/f'{gal}_{key}_maps.npz')
            header=fits.Header.fromstring(str(sample['header']))
            scale=float(np.sqrt(abs(np.linalg.det(WCS(header).celestial.pixel_scale_matrix)))*3600)
            result=segment_variance([sample[f'var_{int(w)}'] for w in WINDOWS],scale,prominence)
            plot_geometry(gal,sample,result,outdir/f'{gal}_geometry.png')
            save_masks(outdir/f'{gal}_approximate_footprints.fits',result,header)
            meta={k:result[k] for k in ['status','n_regions','angle','cuts','local_v_cuts','cells','rectangles','adjacent','pixel_scale','overlap_halfwidth_pixels']}
            meta.update(cube=cube,ao=rec['ao'],windows=list(WINDOWS),prominence=prominence)
            (outdir/f'{gal}_geometry.json').write_text(json.dumps(meta,indent=2)+'\n')
            statuses.append(dict(galaxy=gal,stage='geometry',status=result['status'],regions=result['n_regions']))
            # A detector unable to split a large mosaic must not silently call it single.
            if result['status']!='multiple_candidate':
                print('SKIP RATIOS',gal,result['status'],flush=True);continue
            run_spectra={}
            for rec in rr:
                try:
                    if not rec['continuum']: raise FileNotFoundError('Continuum cube absent')
                    if rec['cube']!=cube: raise ValueError('Run uses a different input cube; geometry needs separate review')
                    if not (float(rec['lmax'])>=float(rec['run'])-100): raise ValueError('CONFIG wavelength range does not match run label')
                    if rec['run']=='8900' and float(rec['lmax'])<8800: raise ValueError('Not an 8900 run')
                    ck=fingerprint(rec['continuum'],json.dumps(meta,sort_keys=True)+str(stride)+'support-v2')
                    sc=cache/f'{gal}_{rec["run"]}_{ck}_spectra.npz'
                    if sc.exists():
                        with np.load(sc) as dd: data={k:dd[k] for k in dd.files}
                    else:
                        data=extract_integrated(rec['continuum'],result,header,rec['ao'],stride)
                        np.savez_compressed(sc,**data)
                    run_spectra[rec['run']]=data
                    statuses.append(dict(galaxy=gal,stage=rec['run'],status='spectra_ok' if np.all(data['counts']>=100) else 'insufficient_support',min_pixels=int(data['counts'].min()),
                                         supported_channels=int(data['supported_channels'].sum()),
                                         min_retained_fraction=float(np.min(data['counts']/np.maximum(1,data['candidate_counts'])))))
                except Exception as exc:
                    statuses.append(dict(galaxy=gal,stage=rec['run'],status='error',error=str(exc)))
                    print('RUN ERROR',gal,rec['run'],repr(exc),flush=True)
            metrics.extend(plot_ratios(gal,result,run_spectra,rr,outdir))
        except Exception as exc:
            statuses.append(dict(galaxy=gal,stage='geometry',status='error',error=str(exc)))
            print('GALAXY ERROR',gal,repr(exc),flush=True)
        finally:
            write_csv(outdir/'status.csv',statuses);write_csv(outdir/'ratio_summary.csv',metrics)
        print('DONE',gal,flush=True)
    return statuses,metrics
