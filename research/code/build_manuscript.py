from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[1];S=ROOT/'manuscript';OUT=ROOT.parent/'output'/'QSEEG'
OUT.mkdir(parents=True,exist_ok=True)
read=lambda n:(S/n).read_text()
parts=[read('preamble_intro.tex'),read('method.tex'),read('theory.tex'),read('simulation_design.tex')]
simulation=read('simulation_results.tex')
simulation=simulation.replace('\\subsection{Sensitivity to ordinary and distribution-specific dependence}',read('table_surface.tex')+'\n'+read('table_counts.tex')+'\n\\subsection{Sensitivity to ordinary and distribution-specific dependence}')
parts+=[simulation,read('realdata_design.tex'),read('table_eeg.tex'),read('realdata_results.tex'),read('conclusion_funding.tex')]
parts+=['\\begingroup\n\\singlespacing\n\\small\n\\setlength{\\bibsep}{7pt}\n'+read('references_build.bbl')+'\n\\endgroup\n\\clearpage\n',read('appendix_proofs.tex')]
app=read('appendix_additional.tex')
for token,file in [('RAW_TABLE','table_raw.tex'),('SMALL_TABLE','table_small.tex'),('STUDENTIZATION_TABLE','table_studentization.tex'),('EEG_SENSITIVITY_TABLE','table_eeg_sensitivity.tex')]:app=app.replace('%% '+token+' %%',read(file))
app=app.replace('[htbp]','[H]')
parts+=[app,'\\end{document}\n']
text='\n'.join(parts)
assert not re.search(r'%% [A-Z_]+ %%',text)
assert not re.search(r'\\(?:input|include|bibliography)\{',text)
assert not re.search(r'(?:Appendix|appendix)[~ ]+\\ref',text)
(OUT/'QSEEG.tex').write_text(text)
intro=text.split('\\section{Introduction}')[1].split('\\section{Method}')[0]
keys={key for group in re.findall(r'\\cite[pt]?\{([^}]+)\}',intro) for key in group.split(',')}
allkeys={key for group in re.findall(r'\\cite[pt]?\{([^}]+)\}',text) for key in group.split(',')}
bibkeys=set(re.findall(r'\\bibitem\[[^\]]*\]\{([^}]+)\}',text))
assert allkeys==bibkeys,(allkeys-bibkeys,bibkeys-allkeys)
print('Introduction references:',len(keys),'total references:',len(bibkeys))
print('Source lines:',len(text.splitlines()))
