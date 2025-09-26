# Slim version of Python
FROM python:3.12

# Switch to new user
USER $USERNAME

# Set working directory and copy app
WORKDIR /APP
COPY APP/ .
COPY requirements.txt .

RUN echo "== /APP ==" && ls -lah /APP && \
    echo "== /APP/INPUTS ==" && ls -lah /APP/INPUTS


# Download Package Information and install dependencies
RUN apt-get update -y && \
    apt-get install -y xauth && \
    apt-get install -y python3-tk tk sudo && \
    rm -rf /var/lib/apt/lists/* && \ 
    pip install --no-cache-dir -r requirements.txt

ARG UID=1000
ARG GID=1000

ENV USERNAME myNewUserName

# Create user, group, sudo permissions
RUN groupadd -g $GID $USERNAME && \
    useradd -m -u $UID -g $GID -s /bin/bash $USERNAME && \
    echo "$USERNAME:$USERNAME" | chpasswd && \
    usermod -aG sudo $USERNAME && \
    mkdir -p /etc/sudoers.d && \
    echo "$USERNAME ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/$USERNAME && \
    chmod 0440 /etc/sudoers.d/$USERNAME

# Run the Tkinter app
CMD ["python3", "Main_inputs_ARTDECO.py"]

