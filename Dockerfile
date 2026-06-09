FROM python:3.12-slim

# Install system dependencies required by RDKit
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    libexpat1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy application source
COPY smarter/ ./smarter/
COPY .streamlit/ ./.streamlit/

EXPOSE 8501

# Streamlit-specific settings: disable the browser auto-open and CORS for server use
ENTRYPOINT ["streamlit", "run", "smarter/app.py", \
    "--server.port=8501", \
    "--server.address=0.0.0.0", \
    "--server.headless=true"]
