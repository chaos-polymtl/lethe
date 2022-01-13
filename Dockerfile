ARG DEALII_IMAGE_VERSION="v9.3.0"

FROM dealii/dealii:${DEALII_IMAGE_VERSION}-focal as builder

# Don't run as dealii user to avoid permission problems
USER root

COPY . lethe

RUN mkdir lethe/build && cd lethe/build && \
    cmake .. \
      -DCMAKE_BUILD_TYPE=Release  \
      -DCMAKE_INSTALL_PREFIX=/tmp/lethe  \
      -DBUILD_SHARED_LIBS=ON && \
    make -j $(nproc) && \
    make install

FROM dealii/dealii:${DEALII_IMAGE_VERSION}-focal as runner

COPY --from=builder /tmp/lethe/bin/ /usr/bin/
COPY --from=builder /tmp/lethe/lib/ /usr/lib/

ENTRYPOINT ["gls_navier_stokes_2d"]