FROM python:3.12-slim AS base

# Install Perl (needed for legacy alignment and Markov scripts)
RUN apt-get update && \
    apt-get install -y --no-install-recommends perl && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy the full repository so that alignment/, markov/, seqtools/ scripts
# are available at runtime.
COPY . .

# Install the Python package
RUN pip install --no-cache-dir .

# Set the repo root so wrappers can locate Perl scripts and scoring
# matrices regardless of how the Python package was installed.
ENV BIOSEA_REPO_ROOT=/app

EXPOSE 8000

CMD ["uvicorn", "api.server:app", "--host", "0.0.0.0", "--port", "8000"]
