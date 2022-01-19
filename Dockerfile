FROM python:3.8

ENV FLASK_APP "adifa"
ENV FLASK_RUN_PORT=5000

RUN mkdir /app
WORKDIR /app

COPY ./requirements.txt /app

RUN pip install --upgrade pip
RUN pip install -r requirements.txt
RUN pip install gunicorn

ADD . /app
ADD https://www.digicert.com/CACerts/BaltimoreCyberTrustRoot.crt.pem /app/BaltimoreCyberTrustRoot.crt.pem

EXPOSE 5000

ENTRYPOINT ["sh", "./boot.sh"]