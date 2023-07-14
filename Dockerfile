ARG DEALII_IMAGE_VERSION="v9.5.0"

FROM dealii/dealii:${DEALII_IMAGE_VERSION}-focal as builder

ARG LETHE_INSTALL_DIR=/opt/lethe

# Don't run as dealii user to avoid permission problems
USER root

# Copy required files and directories for compilation
COPY CMakeLists.txt lethe/
COPY applications_tests lethe/applications_tests
COPY tests lethe/tests
COPY applications lethe/applications
COPY include lethe/include
COPY source lethe/source

# Build
RUN mkdir lethe/build && cd lethe/build && \
    cmake .. \
      -DCMAKE_BUILD_TYPE=Release  \
      -DCMAKE_INSTALL_PREFIX=${LETHE_INSTALL_DIR} \
      -DBUILD_SHARED_LIBS=ON && \
    make -j $(nproc) && \
    make install

FROM dealii/dealii:${DEALII_IMAGE_VERSION}-focal as runner

LABEL org.opencontainers.image.title="lethe" \
      org.opencontainers.image.authors="lethe-cfd" \
      org.opencontainers.image.source="https://github.com/lethe-cfd/lethe"

ARG LETHE_INSTALL_DIR=/opt/lethe

# Set env vars
ENV LETHE_INSTALL_DIR=${LETHE_INSTALL_DIR} \
    PATH=${LETHE_INSTALL_DIR}/bin:${PATH} \
    LD_LIBRARY_PATH=${LETHE_INSTALL_DIRECTORY}/lib

# Copy entrypoint script and make it executable
COPY container/docker_entrypoint.py /usr/local/bin/lethe
RUN sudo chmod +x /usr/local/bin/lethe

# Copy built executables from builder stage
COPY --from=builder ${LETHE_INSTALL_DIR} ${LETHE_INSTALL_DIR}

ENTRYPOINT ["/usr/local/bin/lethe"]
