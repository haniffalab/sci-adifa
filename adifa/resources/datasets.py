from flask import jsonify, request
from flask_restful import Resource

from adifa import models
from adifa.resources.errors import SqlAlchemyError
from adifa.utils import adata_utils


class Datasets(Resource):
    def get(self):
        try:
            return jsonify(datasets=[i.serialize for i in models.Dataset.query.all()])
        except Exception as e:
            raise SqlAlchemyError(e)
        

class Dataset(Resource):
    def get(self, id):
        return models.Dataset.query.get(id).serialize


class Coordinates(Resource):
    def get(self):
        datasetId = request.args.get('datasetId', 0, type=int)
        obsm = request.args.get('embedding', 'X_umap', type=str)

        return adata_utils.get_coordinates(datasetId, obsm)


class Labels(Resource):
    def get(self):
        datasetId = request.args.get('datasetId', 0, type=int)
        obsm = request.args.get('embedding', 'X_umap', type=str)
        gene = request.args.get('gene', '', type=str)
        obs = request.args.get('obs', '', type=str)

        return adata_utils.get_labels(datasetId, obsm, gene=gene, obs=obs)


class Bounds(Resource):
    def get(self):
        datasetId = request.args.get('datasetId', 0, type=int)
        obsm = request.args.get('embedding', 'X_umap', type=str)

        return adata_utils.get_bounds(datasetId, obsm)


class Genelist(Resource):
    def get(self, id):
        searchterm = request.args.get('term', '', type=str)

        return { "searchterm": searchterm, "genes": adata_utils.search_genes(id, searchterm) }

class Genesearch(Resource):
    def get(self, id):
        q = request.args.get('q', '', type=str)

        return { "searchterm": q, "genes": adata_utils.gene_search(id, q) }


class CellByGeneAggregates(Resource):
    def get(self, id):
        genes = request.args.getlist('genes')
        obs = request.args.get('obs')

        output = []

        for gene in genes:
            output += adata_utils.cat_expr_w_counts(id, obs, gene)

        return { "data": output }

class DiseaseGeneList(Resource):
    def get(self, id):
        term = request.args.get('term', '', str)

        output = {}
        disease = 'Haematological abnormality'
        gene = 'Gene mutation'
        datafile = adata_utils.disease_filename()

        #open csv file;
        from csv import DictReader

        with open(datafile, newline='', errors='ignore') as f:
            reader = DictReader(f)
            
            #grab diseases that match search + related genes
            for row in reader:
                if term.lower() in row[disease].lower(): #should also check gene is in dataset here
                    if row[disease] not in output:
                        output[row[disease]] = []
                    
                    output[row[disease]].append(row[gene])
        

        # return resources
        return { 'data': output }