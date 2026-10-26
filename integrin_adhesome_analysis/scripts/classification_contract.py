"""Shared version 1.0 family scope and user-supplied workbook layout contract."""

# Resolve the organized bundle independently of the caller working directory.
import sys as _bundle_sys
from pathlib import Path as _BundlePath
_BUNDLE_ROOT=next(p for p in _BundlePath(__file__).resolve().parents if (p/'pipeline_paths.py').is_file())
if str(_BUNDLE_ROOT) not in _bundle_sys.path:_bundle_sys.path.insert(0,str(_BUNDLE_ROOT))
from pipeline_paths import activate as _activate_bundle_paths
_activate_bundle_paths()
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]


def general_family(name,families):
    aliases={'actin_superfamily_screen':'actin','ADF_H_screen':'cofilin',
        'paxillin_family':'paxillin_like','zyxin_family':'zyxin_like',
        'Src_family':'Src','Csk_family':'Csk','SH3_SH2_TK_screen':'Src'}
    original=name;seen=set()
    while name in families and name not in seen:
        seen.add(name)
        if name=='PINCH_like':return name
        if name in aliases:return aliases[name]
        rule=families[name]
        if rule.get('classification_level')!='subfamily':return name
        parent=rule.get('parent_family')
        if not parent:return name
        name=parent
    return aliases.get(name,original)


def general_rules(config):
    families=config['families']
    selected={name:dict(rule) for name,rule in families.items() if general_family(name,families)==name}
    for name in ('actin','cofilin','paxillin_like','zyxin_like','Src','Csk','PINCH_like'):
        if name in selected:selected[name].update(screening_only=False,classification_level='family')
    return selected


def schema(stage):
    return json.loads((ROOT/'templates/output_schema.json').read_text())['sheets'][str(stage)]


def fields_for(stage,name,rows,reference_species=None,target_species=None):
    spec=next((x for x in schema(stage) if x['name']==name),dict(columns=[]))
    columns=spec['columns']
    if stage==2 and name=='Audited_family_assignment' and reference_species is not None:
        template_species=('Hsap','Cele','Dmel','Mmus','Xlae','Shae','Sjap','Sman')
        dynamic=[]
        for sp in reference_species:
            dynamic.extend('Orthology_'+sp+'_'+suffix for suffix in
                ('direct_orthologues','relationship','HOG_copy_count','HOG_members'))
        for sp in target_species or []:dynamic.append('Orthology_'+sp+'_HOG_copy_count')
        columns=[c for c in columns if not any(c.startswith('Orthology_'+sp+'_') for sp in template_species)]
        at=columns.index('Orthology_reference_support') if 'Orthology_reference_support' in columns else len(columns)
        columns=columns[:at]+dynamic+columns[at:]
    return list(dict.fromkeys([*columns,*(key for row in rows for key in row)]))


def write_sheet(book,stage,name,rows,reference_species=None,target_species=None,fields=None):
    """Production Python export follows the uploaded tabular templates."""
    from openpyxl.styles import Font,PatternFill,Alignment
    from openpyxl.utils import get_column_letter
    sheet=book.create_sheet(name)
    columns=fields or fields_for(stage,name,rows,reference_species,target_species)
    if not columns:columns=['protein_id','evidence_status']
    sheet.append(columns)
    for row in rows:sheet.append([row.get(key) for key in columns])
    spec=next((x for x in schema(stage) if x['name']==name),{})
    for cell in sheet[1]:
        cell.font=Font(name='Arial',size=10,bold=True,color='FFFFFF')
        cell.fill=PatternFill('solid',fgColor='1F4E78');cell.alignment=Alignment(wrap_text=True,vertical='center')
    for i,key in enumerate(columns,1):
        col=get_column_letter(i)
        width=spec.get('widths',{}).get(col,35 if any(x in key.lower() for x in ('architecture','rationale','hypotheses','details')) else 22)
        if stage==2 and name=='Audited_family_assignment':
            width={'Specie':16,'protein_id':28,'putative_audited_family':30,'protein_name':60,'grade':10,'length':12,
                'Pfam_architecture':60,'NCBI_CDD_architecture':80,'architecture_evidence':65,'naming_evidence':65}.get(key,24)
        sheet.column_dimensions[col].width=width
    sheet.row_dimensions[1].height=max(spec.get('header_height') or 45,75 if any(len(key)>60 for key in columns) else 45)
    sheet.freeze_panes='C2' if len(columns)>2 else 'A2'
    sheet.auto_filter.ref=sheet.dimensions;sheet.sheet_view.showGridLines=False
    return sheet


def stage_one_workbook(rows,destination):
    from openpyxl import Workbook
    book=Workbook();book.remove(book.active)
    records=[]
    for row in rows:
        identifier=str(row.get('protein_id',''))
        species=row.get('Specie') or (identifier.split('__',1)[0] if '__' in identifier else '')
        records.append(dict(row,Specie=species))
    write_sheet(book,1,'candidates',records)
    book.save(destination)
