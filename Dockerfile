FROM tiangolo/uwsgi-nginx-flask:python3.8
COPY ./requirements.txt /var/www/requirements.txt
RUN pip install -r /var/www/requirements.txt
COPY ./app/ /app
WORKDIR /app
CMD [ "python", "run.py", "--host=0.0.0.0"]