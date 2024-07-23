FROM openjdk:17-slim
RUN apt-get update && apt-get install -y libxext6
ADD https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar /app/juicer_tools.jar
RUN chmod +x /app/juicer_tools.jar
WORKDIR /app

