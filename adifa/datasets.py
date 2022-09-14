import os
from signal import SIG_DFL

from flask import (
    abort,
    Blueprint,
    current_app,
    flash,
    redirect,
    render_template,
    send_from_directory,
    session,
    request,
    url_for,
)
from sqlalchemy import exc

from adifa import models


bp = Blueprint("datasets", __name__)


def get_modalities(dataset):
    try:
        datasets = models.Dataset.query.where(
            models.Dataset.filename == dataset.filename,
            models.Dataset.modality != dataset.modality,
        ).all()
    except exc.SQLAlchemyError as e:
        abort(500)
    return {d.modality: d.id for d in datasets} if datasets else None


@bp.route("/")
def index():
    return render_template("index.html")


@bp.route("/dataset/<int:id>/scatterplot")
def scatterplot(id):
    try:
        dataset = models.Dataset.query.get(id)
        if not dataset:
            abort(404)
    except exc.SQLAlchemyError as e:
        abort(500)

    authenticated = session.get("auth_dataset_" + str(id), False)
    if dataset.password and not authenticated:
        session["auth_redirect"] = "datasets.scatterplot"
        return redirect(url_for("datasets.password", id=id))

    from collections import OrderedDict
    from operator import getitem

    from collections import defaultdict, OrderedDict
    from operator import getitem

    groups = defaultdict(dict)
    for key, value in dataset.data_obs.items():
        groups[value["group"]][key] = value

    if current_app.config.get("KEEP_OBS_ORDER"):
        pass
        # obs = OrderedDict(dataset.data_obs.items())
    else:
        for key, value in groups.items():
            groups[key] = OrderedDict(
                sorted(value.items(), key=lambda x: getitem(x[1], "name"))
            )

    dataset.other_modalities = get_modalities(dataset)

    return render_template("scatterplot.html", dataset=dataset, obs=groups)


@bp.route("/dataset/<int:id>/matrixplot")
def matrixplot(id):
    try:
        dataset = models.Dataset.query.get(id)
        if not dataset:
            abort(404)
    except exc.SQLAlchemyError as e:
        abort(500)

    # Check protected status
    authenticated = session.get("auth_dataset_" + str(id), False)
    if dataset.password and not authenticated:
        session["auth_redirect"] = "datasets.matrixplot"
        return redirect(url_for("datasets.password", id=id))

    from collections import defaultdict, OrderedDict
    from operator import getitem

    groups = defaultdict(dict)
    for key, value in dataset.data_obs.items():
        groups[value["group"]][key] = value

    if current_app.config.get("KEEP_OBS_ORDER"):
        pass
        # obs = OrderedDict(dataset.data_obs.items())
    else:
        for key, value in groups.items():
            groups[key] = OrderedDict(
                sorted(value.items(), key=lambda x: getitem(x[1], "name"))
            )

    dataset.other_modalities = get_modalities(dataset)

    return render_template("matrixplot.html", dataset=dataset, obs=groups)


@bp.route("/dataset/<int:id>/download", methods=["GET"])
def download(id):
    try:
        dataset = models.Dataset.query.get(id)
        if not dataset:
            abort(404)
    except exc.SQLAlchemyError as e:
        abort(500)

    # Check protected status
    authenticated = session.get("auth_dataset_" + str(id), False)
    if dataset.password and not authenticated:
        session["auth_redirect"] = "datasets.download"
        return redirect(url_for("datasets.password", id=id))

    if dataset.download_link:
        return redirect(dataset.download_link, code=302)
    else:
        return send_from_directory(
            current_app.config.get("DATA_PATH"), dataset.filename, as_attachment=True
        )


@bp.route("/dataset/<int:id>/password", methods=["GET", "POST"])
def password(id):
    try:
        dataset = models.Dataset.query.get(id)
        if not dataset:
            abort(404)
    except exc.SQLAlchemyError as e:
        abort(500)

    authenticated = session.get("auth_dataset_" + str(id), False)
    if dataset.password and authenticated:
        return redirect(
            url_for(session.get("auth_redirect", "datasets.scatterplot"), id=id)
        )

    # Handle the POST request
    if request.method == "POST":
        password = request.form.get("password")
        if dataset.password == password:
            session["auth_dataset_" + str(id)] = True
            return redirect(
                url_for(session.get("auth_redirect", "datasets.scatterplot"), id=id)
            )
        else:
            flash("The password you entered is not correct")

    # Otherwise handle the GET request
    return render_template("password.html", dataset=dataset)
