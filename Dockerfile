# FROM alpine:3.20
FROM ubuntu:20.04

WORKDIR /app
COPY . /app/

RUN apt-get update && apt-get install --no-install-recommends -y \
    g++ \
    valgrind \
    make \
    libgtest-dev \
    lcov \
    gzip && \
    apt-get clean && rm -rf /var/lib/apt/lists/*
EXPOSE 8080
VOLUME [ "/app" ]
CMD ["tail", "-f", "/dev/null"]