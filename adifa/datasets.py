import os

from flask import Blueprint, current_app, send_from_directory, redirect, render_template, session

from adifa import models

bp = Blueprint('datasets', __name__)

@bp.route("/")
def index():
    return render_template('index.html')
  
@bp.route('/dataset/<int:id>/scatterplot')
def scatterplot(id):
    # if session.get("app_password") is None:
    #     return render_template('password.html')
    dataset = models.Dataset.query.get(id)

    from collections import OrderedDict 
    from operator import getitem 
    obs = OrderedDict(sorted(dataset.data_obs.items(), key = lambda x: getitem(x[1], 'name'))) 
    return render_template('scatterplot.html', did=id, dataset=dataset, obs=obs)    
  
@bp.route('/dataset/<int:id>/heatmap')
def heatmap(id):
    dataset = models.Dataset.query.get(id)

    from collections import OrderedDict 
    from operator import getitem 
    obs = OrderedDict(sorted(dataset.data_obs.items(), key = lambda x: getitem(x[1], 'name'))) 
    return render_template('heatmap.html', did=id, dataset=dataset, obs=obs)    

@bp.route('/dataset/<int:id>/download', methods=['GET'])
def download(id):
    dataset = models.Dataset.query.get(id)

    if dataset.download_link:
        return redirect(dataset.download_link, code=302)
    else:
        if (os.path.isabs(current_app.config.get('DATA_PATH'))):
            directory = os.path.realpath(current_app.config.get('DATA_PATH'))
        else:
            directory = os.path.realpath(current_app.root_path + '/../' + current_app.config.get('DATA_PATH'))

        return send_from_directory(
            directory, dataset.filename, as_attachment=True
        )