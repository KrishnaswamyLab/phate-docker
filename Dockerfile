FROM python:3
WORKDIR /app
RUN apt-get update -qq 
RUN pip3 install --upgrade pip
COPY requirements.txt .
RUN pip3 install --prefer-binary -r requirements.txt
COPY run_phate.py .
COPY utils.py .
COPY phate_command_line.py .
ENTRYPOINT [ "python", "./phate_command_line.py" ]
