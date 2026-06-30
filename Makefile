all: rust

# --- Rust build (default) ---
rust:
	cargo build --release --manifest-path evm/Cargo.toml

rust-debug:
	cargo build --manifest-path evm/Cargo.toml

test:
	cargo test --manifest-path evm/Cargo.toml --all
	cd testing && bash runMe.rust.sh

# --- Legacy ParaFly build (only needed for Perl EVM) ---
OS := $(shell uname)

CXX = g++
CC = gcc

parafly:
	cd plugins/ParaFly && sh ./configure --prefix=`pwd` CXX=$(CXX) CC=$(CC) CFLAGS="-fopenmp" CXXFLAGS="-fopenmp" && $(MAKE) install

large_sample_data:
	git clone https://github.com/EVidenceModeler/EVM_sample_data.git
