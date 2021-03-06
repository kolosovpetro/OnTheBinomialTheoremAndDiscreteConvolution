# On the link between binomial theorem and discrete convolution

<p align="center">
  <img src="img/abstract.PNG" alt="logo_example"/>
</p>

[![Build PDF](https://github.com/kolosovpetro/OnTheBinomialTheoremAndDiscreteConvolution/actions/workflows/build-pdf.yml/badge.svg)](https://github.com/kolosovpetro/OnTheBinomialTheoremAndDiscreteConvolution/actions/workflows/build.yml/badge.svg)
[![Build and Deploy PDF](https://github.com/kolosovpetro/OnTheBinomialTheoremAndDiscreteConvolution/actions/workflows/build-and-deploy-pdf.yml/badge.svg)](https://github.com/kolosovpetro/OnTheBinomialTheoremAndDiscreteConvolution/actions/workflows/build-and-deploy.yml/badge.svg)
![contributors count](https://img.shields.io/github/contributors/kolosovpetro/OnTheBinomialTheoremAndDiscreteConvolution)

## What is all about

Source code of the manuscript entitled "On the link between binomial theorem and discrete convolution" along with
Mathematica programs in order to verify results.

## Build and run in Intellij IDEA

- Install `MikTeX`: https://miktex.org/download
- Update `MikTeX`
- Install `SumatraPDF` viewer: https://www.sumatrapdfreader.org/download-free-pdf-viewer
- Install `Intellij IDEA Ultimate`: https://www.jetbrains.com/idea/download/#section=windows
- Activate `Intellij IDEA Ultimate`
- Install `TeXiFy IDEA` plugin: https://plugins.jetbrains.com/plugin/9473-texify-idea
- Clone this repository locally: `https://github.com/kolosovpetro/OnTheBinomialTheoremAndDiscreteConvolution.git`
- Open `github-latex-template` folder in `Intellij IDEA Ultimate` and configure as follows
    - LaTeX Configuration
      ![LaTeX Configuration](img/latex_configuration.PNG?raw=true "LaTeX Configuration")
    - BibTeX Configuration
      ![BibTeX Configuration](img/bibtex_configuration.PNG?raw=true "BibTeX Configuration")
- Configure Inverse Search in `Intellij IDEA` for SumatraPDF: `Tools -> LaTeX -> Configure Inverse Search`
- Compile document using `Shift + F10`

## How to use Mathematica package

- Open the package file `OnTheBinomialTheoremAndDiscreteConvolution.m` in Wolfram Mathematica, I use version 13.0
- Execute the package using `Shift+Enter`
- Open the notebook file `OnTheBinomialTheoremAndDiscreteConvolution.nb`
- Execute the line: `Needs["OnTheBinomialTheoremAndDiscreteConvolution"]`
- Continue your work as desired

## Configure CI / CD

Set repository secrets

- `GH_ACCESS_TOKEN`: Generate Github Personal access token at
  `Settings -> Developer Settings -> Personal access tokens -> Generate mew token` and assign in to
  secret `GH_ACCESS_TOKEN`
- `GH_NAME`: Your Github username
- `GH_EMAIL`: Your Github email

## Actions and their trigger policy

- `build-pdf.yml` builds project using `TeXLive`. Triggered on `pull_request`, `push` to `develop` branch
- `build-and-deploy-pdf.yml` builds project using `TeXLive` and deploys to `GitHub Pages`. Triggered on `push` to `main`
  branch
