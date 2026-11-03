"""Raw-input registry and explicit unidentified_input decisions; standard library only.

Aliases name configured datasets, never alternate historical result directories.
Missing mandatory data always blocks execution; clean-state assessment records
that block without raising so every conditional stage can be inspected.
"""
import json
import hashlib
from pathlib import Path, PureWindowsPath

ROOT = Path(__file__).resolve().parents[2]
PROTEIN_EXTENSIONS = ('.fa', '.faa', '.fasta', '.fas', '.fsa', '.pep', '.aa')


def deployment_source_digest(path):
    """Release-index text normalization only; scientific input hashes stay bytewise."""
    path=Path(path)
    data=path.read_bytes()
    if path.suffix.lower() in ('.py','.md','.yml','.yaml','.txt','.json') or path.name in ('LICENSE','VERSION','.gitattributes','.gitignore'):
        data=data.replace(b'\r\n',b'\n')
    return hashlib.sha256(data).hexdigest()


def is_protein_path(path):
    name = Path(path).name.lower()
    return name.endswith(PROTEIN_EXTENSIONS) or name.endswith(tuple(x + '.gz' for x in PROTEIN_EXTENSIONS))


class InputDependencyError(ValueError):
    def __init__(self, event):
        self.event = event
        super().__init__(f"unidentified_input: stage={event['stage']}, input={event['requested']}, "
                         f"decision={event['decision']}: {event['reason']}. {event['remedy']}")


class InputRegistry:
    def __init__(self, root=ROOT, path=None):
        self.root = Path(root).resolve()
        self.path = Path(path) if path else self.root / 'config/input_registry.json'
        self.config = json.loads(self.path.read_text(encoding='utf-8-sig'))
        if self.config.get('version') != '1.0':
            raise ValueError('Input registry version must be 1.0')
        policy = self.config.get('unidentified_input', {})
        if policy.get('mandatory_missing') != 'stop' or policy.get('optional_missing') not in ('skip', 'stop'):
            raise ValueError('Mandatory input failures must stop; optional_missing must be skip or stop')
        self.aliases = {}
        for key, entry in self.config['datasets'].items():
            raw = entry['path']
            if Path(raw).is_absolute() or PureWindowsPath(raw).is_absolute():
                raise ValueError('Registry paths must be project-root-relative: ' + key)
            if entry['kind'] not in ('file', 'directory', 'protein_file', 'protein_directory'):
                raise ValueError('Unknown input kind: ' + key)
            for name in [key, *entry.get('aliases', [])]:
                if name in self.aliases:
                    raise ValueError('Ambiguous input registry alias: ' + name)
                self.aliases[name] = key
        for name, files in self.config.get('reference_panels', {}).items():
            if not isinstance(files, list) or not files:
                raise ValueError('A reference panel must list explicit nonempty relative FASTA paths: ' + name)
            self.select_references(panel=name, check_exists=False)
        # Keep target streams physically distinct from evaluation inputs.
        target = self.path_for('target_raw')
        reference = self.path_for('reference_raw')
        if target == reference or target in reference.parents or reference in target.parents:
            raise ValueError('Target and reference raw roots must be separate')
        for name in ('calibration_fasta', 'calibration_labels', 'validation_fasta', 'validation_labels'):
            path = self.path_for(name)
            if path == target or target in path.parents:
                raise ValueError('Evaluation controls must be outside target_raw: ' + name)

    def path_for(self, name):
        key = self.aliases.get(name)
        if key is None:
            raise KeyError('Unregistered dataset: ' + name)
        return (self.root / self.config['datasets'][key]['path']).resolve()

    def assess(self, name, stage, mandatory=True, override=None, clean_state=False, execution=False, allow_file=False):
        key = self.aliases.get(name)
        path = Path(override).resolve() if override is not None else self.path_for(name) if key else None
        kind = self.config['datasets'][key]['kind'] if key else 'file'
        reason = 'unregistered_dataset' if key is None and override is None else 'expected_absent_in_clean_state' if clean_state else ''
        available = False
        if not reason:
            if not path.exists():
                reason = 'path_missing'
            elif kind=='protein_directory' and allow_file and path.is_file() and not is_protein_path(path):
                reason = 'unsupported_protein_format'
            elif kind.endswith('directory') and not path.is_dir() and not (allow_file and kind=='protein_directory' and path.is_file() and is_protein_path(path)):
                reason = 'expected_directory'
            elif kind.endswith('file') and not path.is_file():
                reason = 'expected_file'
            elif kind == 'protein_file' and not is_protein_path(path):
                reason = 'unsupported_protein_format'
            elif kind == 'protein_directory' and not (allow_file and path.is_file()) and not any(p.is_file() and is_protein_path(p) for p in path.rglob('*')):
                reason = 'no_supported_protein_fastas'
            else:
                available = True
        decision = 'use' if available else 'stop_stage' if mandatory or self.config['unidentified_input']['optional_missing'] == 'stop' else 'skip_operation'
        event = dict(action='identified_input' if available else 'unidentified_input', requested=name, dataset=key,
                     stage=stage, mandatory=mandatory, path=str(path) if path else None, available=available,
                     decision=decision, reason=reason or 'registered_input_available',
                     remedy='Supply the matching dataset at the registered path, or pass an explicit input path; never use cached results as raw inputs.')
        if key=='target_raw':
            event['remedy']='Supply matching raw target proteomes; use --input-dir <directory> (homology: --target-dir). Never substitute evaluation controls or cached outputs.'
        elif key=='reference_raw':
            event['remedy']='Supply matching raw references; use --reference-raw <directory> (homology: --background-dir). Select only the required reference panel.'
        if execution and decision == 'stop_stage':
            raise InputDependencyError(event)
        return event

    def require(self, name, stage, override=None, allow_file=False):
        event = self.assess(name, stage, override=override, execution=True,allow_file=allow_file)
        return Path(event['path'])

    def select_references(self, panel=None, files=None, raw=None, check_exists=True):
        if panel and files:
            raise ValueError('Choose a registry reference panel or explicit reference files, not both')
        selected = self.config.get('reference_panels', {}).get(panel) if panel else files
        if panel and selected is None:
            raise ValueError('Unidentified reference panel: ' + panel)
        if selected is None:
            return None
        base = Path(raw).resolve() if raw else self.path_for('reference_raw')
        result = []
        for value in selected:
            rel = Path(value)
            if rel.is_absolute() or PureWindowsPath(value).is_absolute() or '..' in rel.parts:
                raise ValueError('Reference selections must be paths within reference_raw: ' + value)
            path = (base / rel).resolve()
            if base not in path.parents or not is_protein_path(path):
                raise ValueError('Invalid reference selection: ' + value)
            if check_exists and not path.is_file():
                event = dict(action='unidentified_input', requested=value, dataset='reference_raw', stage='reference_selection',
                             mandatory=True, path=str(path), available=False, decision='stop_stage', reason='selected_reference_missing',
                             remedy='Supply this selected FASTA or correct the explicit panel; no unrelated reference is substituted.')
                raise InputDependencyError(event)
            result.append(path)
        if not result or len(set(result)) != len(result):
            raise ValueError('Reference selection must be nonempty and unique')
        return result


def raw_input(name):
    return InputRegistry().path_for(name)
