from r-base:4.3.0

# Copy additional scripts from bin and add to PATH
RUN mkdir /opt/bin
COPY *.R .
RUN chmod a+x *.R