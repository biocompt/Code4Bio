FROM shinejh0528/snakemake:7.32.4

# Verificar versión de Python
RUN python --version

# Directorio de trabajo
WORKDIR /app

# Copiar solo lo esencial
COPY Snakefile .
COPY utilities/ ./utilities/

# Configurar permisos
RUN chmod -R a+rwx /app

# Punto de entrada
ENTRYPOINT ["snakemake"]
