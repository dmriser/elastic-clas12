FROM rootproject/root-ubuntu16 

# Run as superuser 
USER root 

#ARG username=user
#RUN userdel builder && useradd --create-home --home-dir /home/${username} ${username}
#ENV home /home/${username}

# Switch to new user 
#USER ${username}

# Setup packages
RUN apt-get update -y
RUN apt-get install -y python-pip
RUN sudo -H pip install numpy==1.16.5 scipy==1.2.2 matplotlib==2.2.4

WORKDIR /data

CMD ["bash", "makePlots.sh"]
# CMD ["python", "mon.py", "-i=rga-out.root", "-o=rga-out"]
# CMD ["python", "mon.py", "-i=rga-in.root", "-o=rga-in"]
# CMD ["python", "mon.py", "-i=sim0.root", "-o=sim0"]
# CMD ["python", "test.py"]