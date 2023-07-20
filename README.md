# Seascapes

This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript on a local computer.

The experiments are meant to run on Linux/Unix/macOS operating systems.

If problems and/or questions are encountered, feel free to [open issues](https://github.com/theoschneider/seascapes/issues).

## 0. Local copy
Clone the repository and `cd` to the dir.
```
git clone https://github.com/theoschneider/seascapes
cd seascapes
```

## 1. Installation

### General dependencies

Install python3 packages.
```
sudo apt install -qq -y python3-dev python3-pip snakemake
pip3 install --user scipy numpy matplotlib pandas statsmodels
```

## 2. Run global analysis


In root folder run `snakemake`:
```
snakemake -k -j 8
```

## 3. Add features
You made modifications to one of the python script, a notebook, this README.md, or you added new features?
You wish this work benefits to all (future) users of this repository?
Please, feel free to open a [pull-request](https://github.com/theoschneider/seascapes/pulls).

## Licence

The MIT License (MIT)

Copyright (c) 2023 Théo Schneider

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


