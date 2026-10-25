"""Execute standard notebook cells in-process; no Jupyter server dependency.

Captures stdout, tables and explicit figure displays into the .ipynb, saves after
every cell (including failures), and propagates execution errors to Slurm.
"""
import sys,json,io,contextlib,base64,traceback,time
from pathlib import Path

path=Path(sys.argv[1]);book=json.loads(path.read_text());namespace={'__name__':'__main__'}
for i,cell in enumerate(book['cells'],1):
    if cell['cell_type']=='code':
        src=cell['source'];src=''.join(src) if isinstance(src,list) else src
        compile(src,f'{path.name}:cell{i}','exec')
def save(): path.write_text(json.dumps(book,indent=1)+'\n')
for cell in book['cells']:
    if cell['cell_type']=='code': cell['outputs']=[];cell['execution_count']=None
save()
count=0
for cell in book['cells']:
    if cell['cell_type']!='code': continue
    count+=1;cell['execution_count']=count;outputs=cell['outputs']
    def png(p): outputs.append(dict(output_type='display_data',data={'image/png':base64.b64encode(Path(p).read_bytes()).decode()},metadata={}))
    def table(frame): outputs.append(dict(output_type='display_data',data={'text/html':frame.to_html(index=False),'text/plain':frame.to_string(index=False)},metadata={}))
    namespace['_record_display_png']=png;namespace['_record_table']=table
    class Tee:
        def __init__(self,name,terminal): self.name=name;self.terminal=terminal
        def write(self,text):
            self.terminal.write(text);self.terminal.flush()
            if outputs and outputs[-1].get('output_type')=='stream' and outputs[-1]['name']==self.name: outputs[-1]['text']+=text
            else: outputs.append(dict(output_type='stream',name=self.name,text=text))
            return len(text)
        def flush(self): self.terminal.flush()
    print('EXECUTING CELL',count,time.strftime('%Y-%m-%d %H:%M:%S'),flush=True)
    try:
        source=cell['source'];source=''.join(source) if isinstance(source,list) else source
        with contextlib.redirect_stdout(Tee('stdout',sys.stdout)),contextlib.redirect_stderr(Tee('stderr',sys.stderr)):
            exec(compile(source,f'{path.name}:cell{count}','exec'),namespace)
    except Exception as exc:
        outputs.append(dict(output_type='error',ename=type(exc).__name__,evalue=str(exc),traceback=traceback.format_exc().splitlines()))
        save();raise
    save()
print('NOTEBOOK EXECUTION COMPLETE',flush=True)
