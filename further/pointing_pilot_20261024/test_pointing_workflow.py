"""Scientific invariants for the reusable notebook helpers."""
import unittest
import numpy as np
import pointing_workflow as pw
import tempfile
from pathlib import Path
from astropy.io import fits

class WorkflowTests(unittest.TestCase):
    def test_ao_window_is_excluded(self):
        wave=np.array([5500.,5770.,5890.,6050.,6400.])
        np.testing.assert_array_equal(pw.allowed_wavelength(wave,True),[True,False,False,False,True])

    def test_shared_overlap_from_six_tiles(self):
        y,x=np.indices((640,440))
        valid=(x>10)&(x<430)&(y>10)&(y<630)
        var=30+3*np.exp(-((x-210)**2+(y-320)**2)/10000)
        for seam in [215]: var=np.where(abs(x-seam)<7,var/2,var)
        for seam in [215,425]: var=np.where(abs(y-seam)<7,var/2,var)
        var=np.where(valid,var,np.nan)
        result=pw.segment_variance([var]*4,pixel_scale=.3)
        self.assertEqual(result['n_regions'],6)
        members=result['footprints'].sum(axis=0)
        self.assertTrue(np.any(members>=2))
        self.assertTrue(np.any(members>=4))

    def test_single_smooth_pointing(self):
        y,x=np.indices((230,230))
        v=20+10*np.exp(-((x-110)**2+(y-110)**2)/1500)
        v=np.where((x>10)&(x<220)&(y>10)&(y<220),v,np.nan)
        self.assertEqual(pw.segment_variance([v]*4,pixel_scale=.3)['n_regions'],1)

    def test_staggered_column_refines_boundary(self):
        y,x=np.indices((650,450))
        valid=(x>10)&(x<440)&(y>10)&(y<640)
        v=np.full(x.shape,30.)
        v[abs(x-225)<7]/=2
        v[(x<225)&(abs(y-320)<7)]/=2
        v[(x>=225)&(abs(y-375)<7)]/=2
        v=np.where(valid,v,np.nan)
        r=pw.segment_variance([v]*4,pixel_scale=.3)
        self.assertTrue(any(abs(c-375)<10 for c in r.get('local_v_cuts',{}).get(1,[])))

    def test_normalized_ratio_and_invalid_denominator(self):
        w=np.array([5100,5400,5500,5600,6500.])
        a=np.array([2,4,6,8,10.]); b=a*3
        r=pw.normalized_ratio(w,a,b,norm_band=(5300,5700))
        np.testing.assert_allclose(r,1)
        b[0]=0
        self.assertTrue(np.isnan(pw.normalized_ratio(w,a,b,norm_band=(5300,5700))[0]))

    def test_fixed_spatial_support_ignores_only_ao_gap(self):
        h=fits.Header()
        for key,val in dict(CTYPE1='RA---TAN',CTYPE2='DEC--TAN',CTYPE3='AWAV',CUNIT1='deg',CUNIT2='deg',CUNIT3='Angstrom',CRPIX1=1.,CRPIX2=1.,CRPIX3=1.,CRVAL1=180.,CRVAL2=10.,CRVAL3=5500.,CDELT1=-.2/3600,CDELT2=.2/3600,CDELT3=200.).items():h[key]=val
        cube=np.ones((5,12,12),dtype='f4')
        cube[1,0,0]=np.nan # 5700: reject this spaxel at EVERY output channel.
        cube[2,:,:]=np.nan # 5900 AO gap: do not reject the other 143 spaxels.
        mask=np.ones((12,12),bool)
        result=dict(labels=mask.astype(int),interiors=[mask])
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'cube.fits';fits.PrimaryHDU(cube,header=h).writeto(p)
            data=pw.extract_integrated(p,result,h,ao=True,stride=1)
        self.assertEqual(int(data['counts'][0]),143)
        self.assertTrue(np.isnan(data['spectra'][0,2]))
        np.testing.assert_allclose(data['spectra'][0,[0,1,3,4]],143)

    def test_globally_unsupported_endpoint_does_not_empty_regions(self):
        h=fits.Header()
        for key,val in dict(CTYPE1='RA---TAN',CTYPE2='DEC--TAN',CTYPE3='AWAV',CUNIT1='deg',CUNIT2='deg',CUNIT3='Angstrom',CRPIX1=1.,CRPIX2=1.,CRPIX3=1.,CRVAL1=180.,CRVAL2=10.,CRVAL3=5500.,CDELT1=-.2/3600,CDELT2=.2/3600,CDELT3=200.).items():h[key]=val
        cube=np.ones((5,12,12),dtype='f4');cube[-1]=np.nan
        mask=np.ones((12,12),bool);result=dict(labels=mask.astype(int),interiors=[mask])
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'cube.fits';fits.PrimaryHDU(cube,header=h).writeto(p)
            data=pw.extract_integrated(p,result,h,stride=1)
        self.assertEqual(int(data['counts'][0]),144)
        np.testing.assert_allclose(data['spectra'][0,:-1],144)
        self.assertTrue(np.isnan(data['spectra'][0,-1]))

if __name__=='__main__': unittest.main()
