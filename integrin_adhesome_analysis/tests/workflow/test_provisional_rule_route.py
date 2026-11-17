"""Synthetic routing tests, with no ML, BLAST or biological analysis dependencies."""
import ast
from pathlib import Path
import unittest

SOURCE=Path(__file__).resolve().parents[2]/'resource_library/adhesome_ml/family_gate.py'
TREE=ast.parse(SOURCE.read_text())
NAMES={'provisional_rule_resolution','rule_route_exclusions','rule_evaluation','family_identity_rule_proposals'}
def ancestors(name, parents):
    result=[]
    while name in parents:
        name=parents[name]; result.append(name)
    return result
def resolve(level,family,reason,parents,policy):
    return (level,family,reason,'core',family) if family in policy['allowed'] else None
NS={'ancestors':ancestors,'_resolve_family_call':resolve,'canonical':lambda x,a:a.get(x,x)}
exec(compile(ast.Module(body=[n for n in TREE.body if isinstance(n,ast.FunctionDef) and n.name in NAMES],type_ignores=[]),str(SOURCE),'exec'),NS)

class RoutingTests(unittest.TestCase):
    def route(self,families,ranking=(('A',.1),),membership='adhesome_associated',raw=.8,excluded=()):
        return NS['provisional_rule_resolution']([{'family':x} for x in families],membership,raw,.5,ranking,{'A':'parent','B':'parent'},{'allowed':['A','B','parent']},excluded)
    def test_concordant_low_vote(self): self.assertEqual(self.route(['A'])[1:3],('A','provisional_rule_rf_concordant'))
    def test_discordant_low_vote(self): self.assertEqual(self.route(['B'])[1:3],('B','provisional_rule_rf_discordant'))
    def test_duplicate_same_family(self): self.assertEqual(self.route(['A','A'])[1],'A')
    def test_conflicting_shared_parent(self): self.assertEqual(self.route(['A','B'])[:2],('broad','parent'))
    def test_unrelated_conflicts(self): self.assertIsNone(self.route(['A','C']))
    def test_top_three_alone_cannot_call(self): self.assertIsNone(self.route([],ranking=(('A',.4),('B',.3),('C',.3))))
    def test_membership_and_exclusion(self):
        self.assertIsNone(self.route(['A'],membership='non_adhesome'))
        self.assertIsNone(self.route(['A'],raw=.1))
        self.assertIsNone(self.route(['A'],excluded=['background']))
    def test_explicit_exclusion_and_no_label_leakage(self):
        policy={'sequence_exclusions':[{'id':'background','all_features':{'cdd:bad':1}}]}
        self.assertEqual(NS['rule_route_exclusions']({'cdd:bad':1},policy),['background'])
        self.assertEqual(NS['rule_route_exclusions']({'expected_label':'negative'},policy),[])
    def test_insufficient_and_curated_evidence(self):
        rule={'id':'A','family':'A','required_groups':[{}],'specific_cdd_pattern':'A'}
        vector={'identity_rule:A:required_all_satisfied':1}
        self.assertEqual(NS['family_identity_rule_proposals']({},vector,[rule]),[])
        vector['identity_rule:A:specific_cdd_pattern_present']=1
        self.assertEqual(len(NS['family_identity_rule_proposals']({},vector,[rule])),1)
        rule['require_ncbi_curated_specific_hit']=True
        self.assertEqual(NS['family_identity_rule_proposals']({},vector,[rule]),[])
if __name__=='__main__': unittest.main()
