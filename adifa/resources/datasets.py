from collections import OrderedDict

from flask import jsonify, request
from flask_restful import Resource

from adifa import models
from adifa.utils import adata_utils


class Datasets(Resource):
    def get(self):
        return jsonify(datasets=[i.serialize for i in models.Dataset.query.all()])


class Dataset(Resource):
    def get(self, id):
        return models.Dataset.query.get(id).serialize


class Coordinates(Resource):
    def get(self):
        datasetId = request.args.get("datasetId", 0, type=int)
        obsm = request.args.get("embedding", "X_umap", type=str)

        return adata_utils.get_coordinates(datasetId, obsm)


class Masks(Resource):
    def get(self):
        datasetId = request.args.get("datasetId", 0, type=int)

        return adata_utils.get_masks(datasetId)


class Labels(Resource):
    def get(self):
        datasetId = request.args.get("datasetId", 0, type=int)
        obsm = request.args.get("embedding", "X_umap", type=str)
        gene = request.args.get("gene", "", type=str)
        obs = request.args.get("obs", "", type=str)

        return adata_utils.get_labels(datasetId, obsm, gene=gene, obs=obs)


class Bounds(Resource):
    def get(self):
        datasetId = request.args.get("datasetId", 0, type=int)
        obsm = request.args.get("embedding", "X_umap", type=str)

        return adata_utils.get_bounds(datasetId, obsm)


class SearchGenes(Resource):
    def get(self, id):
        q = request.args.get("search", "", type=str)

        output = []
        for gene in adata_utils.search_genes(id, q):
            sample = {"id": gene, "text": gene}
            output.append(sample)

        return {"results": output[:30]}


class SearchDiseases(Resource):
    def get(self, id):
        q = request.args.get("search", "", type=str)

        output = {}
        categories = {}
        disease = "Disease"
        gene = "Gene mutation"
        category = "Category"
        info = "Info"
        datafile = adata_utils.disease_filename()

        # open csv file;
        from csv import DictReader

        with open(datafile, newline="", errors="ignore") as f:
            reader = DictReader(f)

            # grab diseases that match search + related genes
            for row in reader:
                if (
                    q.lower() in row[disease].lower()
                ):  # should also check gene is in dataset here
                    if row[disease] not in output:
                        output[row[disease]] = {}
                        categories[row[disease]] = set()

                    cat = row.get(category, "default").replace("_", " ").upper()
                    categories[row[disease]].add(cat)
                    output[row[disease]].setdefault(cat, {}).setdefault(
                        row[gene], {"gene": row[gene], "info": []}
                    )
                    if info in row:
                        output[row[disease]][cat][row[gene]]["info"].append(row[info])

        # return output
        results = []
        for key, value in output.items():
            sample = {
                "id": key,
                "text": key,
                "values": {
                    k: OrderedDict(sorted(v.items())) for k, v in sorted(value.items())
                },
                "categories": list(categories[key]),
            }
            results.append(sample)

        return {"results": results}


class CellByGeneAggregates(Resource):
    def get(self, id):
        genes = request.args.getlist("genes")
        obs = request.args.get("obs")

        output = []

        for gene in genes:
            output += adata_utils.cat_expr_w_counts(id, obs, gene)

        return {"data": output}
