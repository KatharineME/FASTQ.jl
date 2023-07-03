## Compile fastp from source

```bash
brew install isa-l

brew install libdeflate

git clone https://github.com/OpenGene/fastp.git

cd fastp

make

sudo make install

```

## Compile minimap2 from source

```bash
# Download source code for minimap2-2.24 from https://github.com/lh3/minimap2/releases

cd minimap2-2.24

make
```

## Compile star from source

```bash
# Download STAR 2.7.9.a from https://github.com/alexdobin/STAR/releases

brew install gcc

cd STAR-2.7.9a/source

make STARforMacStatic CXX=/usr/local/Cellar/gcc/11.2.0_3/bin/g++-11

# STAR executable will be in STAR-2.7.9a/bin/MacOSX_x86_64
```
