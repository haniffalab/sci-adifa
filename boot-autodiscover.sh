#!/bin/bash
flask init-db
flask autodiscover
exec gunicorn -b :5000 --timeout=900 --limit-request-line=8190 --access-logfile - --error-logfile - startup:app