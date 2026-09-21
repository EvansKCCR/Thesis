import json

class AdhesomeLibrary:
    def __init__(self, json_filepath):
        """Initialize the library by loading the JSON schema."""
        with open(json_filepath, 'r') as file:
            self.data = json.load(file)['integrin_adhesome_network']

    def get_categories(self):
        """Return a list of all biological categories in the network."""
        return [category['category'] for category in self.data]

    def search_by_protein(self, protein_name):
        """Retrieve all motifs and domains associated with a specific protein."""
        results = []
        for category in self.data:
            for entity in category['entities']:
                if protein_name.lower() in entity['name'].lower():
                    results.append({
                        "Protein": entity['name'],
                        "Category": category['category'],
                        "Modules": entity['modules']
                    })
        return results

    def find_module_by_type(self, module_type):
        """Filter the network by module type (e.g., 'Motif' or 'Domain')."""
        results = []
        for category in self.data:
            for entity in category['entities']:
                matched_modules = [m for m in entity['modules'] if m['type'].lower() == module_type.lower()]
                if matched_modules:
                    results.append({
                        "Protein": entity['name'],
                        "Modules": matched_modules
                    })
        return results

# Example Usage:
# library = AdhesomeLibrary('adhesome_motifs.json')
# print(library.search_by_protein('Talin'))
# print(library.find_module_by_type('Motif'))