#!/bin/bash
flask init-db
exec gunicorn -b :5000 --timeout=900 --limit-request-line=8190 --access-logfile - --error-logfile - startup:app