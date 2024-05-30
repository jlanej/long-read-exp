FROM python:3.9-alpine
# install dependencies for the project

RUN apk update && apk upgrade && apk add --no-cache make

RUN pip3 install --upgrade pip
#RUN pip3 install pysam
# install the project dependencies in requirements.txt
COPY requirements.txt /app/requirements.txt
WORKDIR /app
RUN pip3 install -r requirements.txt
# copy the project files into the container
COPY . /app
WORKDIR /app

# run the project
CMD ["python3", "main.py"]