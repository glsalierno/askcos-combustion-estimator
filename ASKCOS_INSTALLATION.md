# ASKCOS Installation Guide for Linux

This guide provides step-by-step instructions to install and run the ASKCOS API locally on a Linux system. ASKCOS (Automated System for Knowledge-based Continuous Organic Synthesis) is an open-source tool developed by MIT for chemical synthesis prediction. Our combustion estimator tool relies on a local ASKCOS server for forward synthesis predictions.

**Important Notes**:
- This setup uses Docker, as recommended by the ASKCOS team for ease and reproducibility.
- Hardware Requirements: At least 16GB RAM, a modern CPU, and sufficient disk space (models can be large). GPU acceleration is optional but recommended for faster predictions.
- Time Estimate: 1-2 hours for initial setup, depending on downloads.
- For official docs, always refer to the **[ASKCOS Documentation](https://askcos.mit.edu/docs/)** and **[ASKCOS GitHub](https://github.com/ASKCOS/askcos-core)**.

## Prerequisites
- **Linux OS**: Tested on Ubuntu 20.04+ or similar (e.g., Debian, Fedora). Use a server-grade setup if running heavy workloads.
- **Docker and Docker Compose**: Install if not already present.
  - Install Docker: Follow **[Docker's Linux Installation Guide](https://docs.docker.com/engine/install/ubuntu/)** (or your distro's equivalent).
  - Install Docker Compose: `sudo apt install docker-compose` (or via pip: `pip install docker-compose`).
- **Git**: `sudo apt install git`.
- **Python 3.8+**: For any local scripting, but not strictly needed for ASKCOS itself.
- **Internet Access**: Required for downloading Docker images and models.

## Step 1: Clone the ASKCOS Repository
Clone the core ASKCOS repository:

git clone https://github.com/ASKCOS/askcos-core.git
cd askcos-core

This contains the Docker Compose files and configuration.

# Example .env
ASKCOS_DB_HOST=askcos-db
ASKCOS_DB_PORT=27017
ASKCOS_DB_NAME=askcos
CELERY_BROKER_URL=redis://askcos-redis:6379/0
CELERY_RESULT_BACKEND=redis://askcos-redis:6379/0

- Adjust ports if there's conflicts (e.g., if 27017 is in use).
- For forward synthesis (used in our tool), ensure the relevant models are enabled in `docker-compose.yml`.

## Step 3: Build and Start the Docker Containers
Run Docker Compose to build and start the services:

docker-compose up -d --build

- `-d`: Runs in detached mode (background).
- `--build`: Builds images if needed.
- This starts services like MongoDB (database), Redis (task queue), Celery workers, and the API server.
- Wait 5-10 minutes for initial startup and model downloads. Check logs with `docker-compose logs -f`.

Key services:
- **askcos-api**: The main API server (accessible at `http://0.0.0.0:8000` or your configured port).
- **askcos-db**: MongoDB for data storage.
- **askcos-celery**: Workers for prediction tasks (e.g., forward synthesis).

## Step 4: Verify the Installation
- Check if containers are running: `docker-compose ps` (all should be "Up").
- Test the API: Open a browser or use curl to hit `http://0.0.0.0:8000/api/health/` (should return a status message).
- For forward synthesis (relevant to our tool): Send a test request via curl:

curl -X POST "http://0.0.0.0:8000/api/forward/call-sync" 
-H "Content-Type: application/json" 
-d '{"smiles": ["CCO"], "reagents": "O=O"}'

- Expect a JSON response with predicted products.

If the API isn't responding:
- Check logs: `docker-compose logs askcos-api`.
- Ensure ports are open (e.g., firewall rules: `sudo ufw allow 8000`).

## Step 5: Integrate with Our Tool
- Once ASKCOS is running, update the base URL in `process_reactivity_oxygen.py` if needed (default: `http://0.0.0.0`).
- Run our scripts as described in the main README.md.

## Troubleshooting
- **Docker Errors**: Ensure you're in the sudo group for Docker (`sudo usermod -aG docker $USER`) and restart your session.
- **Model Download Failures**: If models fail to load, check internet connectivity or increase Docker resource limits (via Docker Desktop settings or daemon.json).
- **Port Conflicts**: Edit `docker-compose.yml` to change exposed ports (e.g., API from 8000 to 8080).
- **Performance Issues**: For faster predictions, enable GPU: Add NVIDIA runtime to Docker (requires CUDA drivers; see **[NVIDIA Docker Guide](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)**).
- **Updates**: Pull latest changes: `git pull` and rebuild with `docker-compose up -d --build`.
- Common Issues: If Celery tasks fail, restart workers: `docker-compose restart askcos-celery`.

## Advanced Configuration
- **Custom Models**: Download additional models from ASKCOS releases and mount them via Docker volumes.
- **Scaling**: For production, scale Celery workers: `docker-compose up -d --scale askcos-celery=4`.
- **Security**: For local use only; do not expose publicly without authentication.

If you encounter issues, check the **[ASKCOS Issues Page](https://github.com/ASKCOS/askcos-core/issues)** or forums.

*Last updated: February 2026*

## Step 2: Configure Environment
- Edit the `.env` file (create one if it doesn't exist) in the root directory. Key settings:
