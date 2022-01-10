import hashlib
import re

from flask import current_app, Response, request, jsonify
from math import isnan
import numpy as np
import scanpy as sc
from scipy import stats
from scipy.sparse import find

from adifa import models

def get_annotations(adata):
	annotations = {'obs': {}, 'obsm': {}}
	
	switcher = {
		'category': type_category,
		'bool': type_bool,
		'int': type_numeric,
		'float': type_numeric,
		'complex': type_numeric,
	}

	for name in adata.obs:
		# Map numpy dtype to a simple type for switching
		dtype = re.sub(r'[^a-zA-Z]', '', adata.obs[name].dtype.name)
		# Get the function from switcher dictionary
		func = switcher.get(dtype, type_discrete)
		# Define an API key safe 
		slug = re.sub(r'[^a-zA-Z0-9]', '', name).lower()
		annotations['obs'][slug] = func(adata.obs[name])
		annotations['obs'][slug]['name'] = name
		
	annotations['obsm'] = [ value for value in adata.obsm ]

	return annotations

def get_bounds(datasetId, obsm):
	dataset = models.Dataset.query.get(datasetId)
	adata = current_app.adata[dataset.filename]
	#adata = current_app.adata
	#adata = sc.read(current_app.config.get('DATA_PATH') + 'covid_portal.h5ad')      

	# Embedded coordinate bounds
	output = {
		'x': {
			'min': adata.obsm[obsm][:,0].min().item(),
			'max': adata.obsm[obsm][:,0].max().item()
		},
		'y': {
			'min': adata.obsm[obsm][:,1].min().item(),
			'max': adata.obsm[obsm][:,1].max().item()
		}
	}

	return output

def get_coordinates(datasetId, obsm):
	dataset = models.Dataset.query.get(datasetId)
	adata = current_app.adata[dataset.filename]
	#adata = current_app.adata
	#adata = sc.read(current_app.config.get('DATA_PATH') + 'covid_portal.h5ad')      

	# True resolution sample generation
	output = []
	for x in adata.obsm[obsm]:
		#output.append(x[:2].tolist())
		output.append([round(num, 4) for num in x[:2].tolist()])

	return output

def get_labels(datasetId, obsm, gene="", obs=""):
	dataset = models.Dataset.query.get(datasetId)
	adata = current_app.adata[dataset.filename]
	#adata = current_app.adata
	#adata = sc.read(current_app.config.get('DATA_PATH') + 'covid_portal.h5ad')      

	if (gene):
		try:
			output = [0] * len(adata.obsm[obsm])
			#expression = adata[:,gene].X/max(1,adata[:,gene].X.max())
			expression = adata[:,gene].X
			(x,y,v) = find(expression)
			for index, i in enumerate(x):
				output[i] = str(round(v[index], 4))	
		except KeyError:
			# @todo HANDLE ERROR
			output = [0] * len(adata.obsm[obsm])
		except IndexError:
			# @todo HANDLE ERROR
			output = [0] * len(adata.obsm[obsm])
	elif (obs):
		output = []
		for index, x in enumerate(adata.obsm[obsm]):
			try:	
				output.append(str(adata.obs[obs][index]))
			except KeyError:
				# @todo HANDLE ERROR
				output.append(0)
			except IndexError:
				# @todo HANDLE ERROR
				output.append(0)

	return output
		
	# Define obsm cache
	cache_key_obsm = 'coords-obsm' + str(datasetId) + obsm
	cache_hash_obsm = hashlib.sha1(cache_key_obsm.encode()).hexdigest()

	# True resolution sample generation
	output['samples'] = []

	# Define obs data lists
	obs_data = {}
	for ob in obs:
		cache_string = 'sample-' + str(datasetId) + '-' + obsm + '-' + ob
		hash_object = hashlib.sha1(cache_string.encode())
		hex_dig = hash_object.hexdigest()
		obs_data[ob] = cache.get(hex_dig)
		if not obs_data[ob]:
			slug = re.sub(r'[^a-zA-Z0-9]', '', ob).lower()
			obs_data[ob] = {'data': [], 'bin_value': []}

			for index, x in enumerate(adata.obsm[obsm]):
				try:
					d = adata.obs[ob][index]
					obs_data[ob]['data'].append(str(d))
					if dataset.data_obs[slug]['type'] == 'categorical':
						# TODO refactor
						if str(adata.obs[ob][index]) in dataset.data_obs[slug]['values'].values():
							v = list(dataset.data_obs[slug]['values'].keys())[list(dataset.data_obs[slug]['values'].values()).index(str(adata.obs[ob][index]))]
						else:
							v = -1
						obs_data[ob]['bin_value'].append(v)
					elif dataset.data_obs[slug]['type'] == 'continuous':
						obs_data[ob]['bin_value'].append(d)	
				except KeyError:
					# TODO handle error
					obs_data[ob]['data'].append(0)
					obs_data[ob]['bin_value'].append(0)
				except IndexError:
					# TODO handle error
					obs_data[ob]['data'].append(0)
					obs_data[ob]['bin_value'].append(0)

			cache.set(hex_dig, obs_data[ob], timeout=0)

	# Define gene data lists
	genes_data = {}
	for gene in genes:
		cache_string = 'sample-gene-' + str(datasetId) + '-' + obsm + '-' + gene
		hash_object = hashlib.sha1(cache_string.encode())
		hex_dig = hash_object.hexdigest()
		genes_data[gene] = cache.get(hex_dig)
		if not genes_data[gene]:
			try:
				adata
			except NameError:
				adata = current_app.adata
			genes_data[gene] = []
			genes_data[gene] = [0] * len(adata.obsm[obsm])
			expression = adata[:,gene].X/max(1,adata[:,gene].X.max())
			(x,y,v) = find(expression)
			for index, i in enumerate(x):
				genes_data[gene][i] = v[index]

			cache.set(hex_dig, genes_data[gene], timeout=0)

	# Append obs and gene data to samples
	for index, x in enumerate(output['samples']):
		for ob in obs:
			output['samples'][index][ob] = obs_data[ob]['data'][index]
		for gene in genes:
			output['samples'][index][gene] = genes_data[gene][index]

	return output

	# Setup bins
	adata_obsm = cache.get(cache_hash_obsm)
	if adata_obsm is None:
		try:
			adata
		except NameError:
			adata = current_app.adata

		adata_obsm = adata.obsm[obsm]

	x = adata_obsm[:,0]
	y = adata_obsm[:,1]
	binx = np.linspace(x.min().item(), x.max().item(), 170)
	biny = np.linspace(y.min().item(), y.max().item(), 170)

	cache_string = 'binned-' + str(datasetId)
	hash_object = hashlib.sha1(cache_string.encode())
	hex_dig = hash_object.hexdigest()
	binned = cache.get(hex_dig)
	if not binned:
		binned = stats.binned_statistic_2d(x, y, None, 'count', bins=[binx, biny])

	binned_stats = {}
	# Get binned statistics for obs
	for ob in obs:
		cache_string = 'binned-' + str(datasetId) + ob
		hash_object = hashlib.sha1(cache_string.encode())
		hex_dig = hash_object.hexdigest()
		binned_stats[ob] = cache.get(hex_dig)
		if not binned_stats[ob]:
			slug = re.sub(r'[^a-zA-Z0-9]', '', ob).lower()
			if dataset.data_obs[slug]['type'] == 'categorical':
				binned_stats[ob] = stats.binned_statistic_2d(x, y, values=obs_data[ob]['bin_value'], statistic=lambda d: mode(d), bins=[binx, biny])
			elif dataset.data_obs[slug]['type'] == 'continuous':
				binned_stats[ob] = stats.binned_statistic_2d(x, y, values=obs_data[ob]['bin_value'], statistic='median', bins=[binx, biny])
			elif dataset.data_obs[slug]['type'] == 'discrete':
				binned_stats[ob] = stats.binned_statistic_2d(x, y, values=obs_data[ob]['bin_value'], statistic='count', bins=[binx, biny])
			
			cache.set(hex_dig, binned_stats[ob], timeout=0)

	# Get binned statistics for genes
	for gene in genes:
		cache_string = 'binned-' + str(datasetId) + gene
		hash_object = hashlib.sha1(cache_string.encode())
		hex_dig = hash_object.hexdigest()
		binned_stats[gene] = cache.get(hex_dig)
		if not binned_stats[gene]:		
			binned_stats[gene] = stats.binned_statistic_2d(x, y, values=genes_data[gene], statistic='median', bins=[binx, biny])

	# Recalculate samples based on binned stats
	rows = binned.statistic.shape[0]
	cols = binned.statistic.shape[1]
	count = 0
	output['samples'].clear()
	for row in range(0, rows):
		for col in range(0, cols):
			insert = False
			sample = {
				"x": binned.x_edge[row],
				"y": binned.y_edge[col],
				"i": count,
			}

			# Insert if bin contains points
			if binned.statistic[row,col]:
				insert = True

			for ob in obs:
				slug = re.sub(r'[^a-zA-Z0-9]', '', ob).lower()
				if binned_stats[ob].statistic[row,col] > 0:
					try:
						if dataset.data_obs[slug]['type'] == 'categorical':
							# TODO review type conversions
							obkey = str(int(binned_stats[ob].statistic[row,col]))
							sample[ob] = dataset.data_obs[slug]['values'].get(obkey, "missing")
						elif dataset.data_obs[slug]['type'] == 'continuous':
							sample[ob] = binned_stats[ob].statistic[row,col]

						insert = True
					except KeyError:
						v = 0
						# TODO handle error

			for gene in genes:
				
				if not isnan(binned_stats[gene].statistic[row,col]):
					try:
						sample[gene] = binned_stats[gene].statistic[row,col]
						insert = True
					except KeyError:
						v = 0
						# TODO HANDLE ERROR						
			if insert:
				output['samples'].append(sample)
				count += 1

	return output['samples']

def search_genes(datasetId, searchterm):
	dataset = models.Dataset.query.get(datasetId)
	adata = current_app.adata[dataset.filename]
	#adata = current_app.adata
	output = [g for g in adata.var_names if searchterm.lower() in g.lower()]

	return output

def gene_search(datasetId, searchterm):
	dataset = models.Dataset.query.get(datasetId)
	adata = current_app.adata[dataset.filename]
	#adata = current_app.adata
	genes = [g for g in adata.var_names if searchterm in g]

	output = []
	for gene in genes:
		sample = {
			"name": gene
		}
		output.append(sample)

	return output

def categorised_expr(datasetId, cat, gene, func="mean"):
	dataset = models.Dataset.query.get(datasetId)
	adata = current_app.adata[dataset.filename]

	data = adata[:,[gene]].to_df()
	grouping = data.join(adata.obs[cat]).groupby(cat)

	if func == "mean":
		expr = grouping.mean()
	elif func == "median":
		expr = grouping.median()

	#counts = grouping.count()/grouping.count().sum()
	#'count': counts.loc[group,gene]
	output = [ {'gene': gene, 'cat': group, 'expr': float(expr.loc[group,gene])} for group in grouping.groups.keys()]

	return output

def cat_expr_w_counts(datasetId, cat, gene, func="mean"):
	from numpy import NaN

	dataset = models.Dataset.query.get(datasetId)
	adata = current_app.adata[dataset.filename]

	groupall = adata[:,[gene]].to_df().join(adata.obs[cat]).groupby(cat)
	groupexpr = adata[:,[gene]].to_df().replace(float(adata[:,[gene]].X.min()), NaN).join(adata.obs[cat]).groupby(cat)

	if func == "mean":
		expr = groupexpr.mean()
	elif func == "median":
		expr = groupexpr.median()

	countpc = (groupexpr.count()*100/groupall.count()).astype(int)

	output = [ {'gene': gene, 'cat': group, 'expr': float(expr.loc[group,gene]), 'count': int(countpc.loc[group,gene])} for group in groupall.groups.keys()]

	return output


def mode(d):
	from scipy import stats as s
	mode = s.mode(d)
	return int(mode[0])

def type_category(obs):
	return { 'type': 'categorical', 'values': dict(enumerate(obs.cat.categories.values.flatten(), 1)) }

def type_bool(obs):
	return { 'type': 'categorical', 'values': ['True','False']}

def type_numeric(obs):
	accuracy = 4
	return { 
		'type': 'continuous',
		'min': round(series_min(obs), accuracy),
		'max': round(series_max(obs), accuracy),
		'mean': round(series_mean(obs), accuracy),
		'median': round(series_median(obs), accuracy)
	}

def type_discrete(obs):
	return { 'type': 'discrete' }

def series_max(s):
	return 0
	if s.isna().all():
		return 0
	else:
		return s.max().item()

def series_min(s):
	return 0
	if s.isna().all():
		return 0
	else:
		return s.min().item()

def series_mean(s):
	return 0
	if s.isna().all():
		return 0
	else:
		return s.mean().item()

def series_median(s):
	return 0
	if s.isna().all():
		return 0
	else:
		return s.median().item()

def disease_filename():
	return current_app.config.get('DATA_PATH') + 'fbm_disease_data.csv'