# Static qqman binary in an empty image.
#
#   docker build -t qqman .
#   docker run --rm qqman                # prints qqman --help
#   docker run --rm -v "$PWD:/data" -w /data qqman qqman manhattan results.tsv -o manhattan.png
FROM rust:1.94-alpine AS build
RUN apk add --no-cache musl-dev
WORKDIR /src/qqman
COPY Cargo.toml Cargo.lock NOTICE ./
COPY src ./src
RUN cargo build --release --locked \
 && strip target/release/qqman

FROM scratch
COPY --from=build /src/qqman/target/release/qqman /usr/local/bin/qqman
COPY LICENSE NOTICE /usr/share/doc/qqman/
CMD ["qqman", "--help"]
