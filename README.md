# On the link between binomial theorem and discrete convolution

Let $\mathbf{P}^{m}_{b}(x)$ be a $2m+1$-degree polynomial in $x$ and $b \in \mathbb{R}$,
\[
\mathbf{P}^{m}_{b}(x) = \sum_{k=0}^{b-1} \sum_{r=0}^{m} \mathbf{A}_{m,r} k^r (x-k)^r,
\]
where $\mathbf{A}_{m,r}$ are real coefficients.
In this manuscript, we introduce the polynomial $\mathbf{P}^{m}_{b}(x)$ and study its properties,
establishing a polynomial identity for odd-powers in terms of this polynomial.
Based on mentioned polynomial identity for odd-powers,
we explore the connection between the Binomial theorem and discrete convolution of odd-powers,
further extending this relation to the multinomial case.
All findings are verified using Mathematica programs.

## Build and run in Intellij IDEA

- Install `MikTeX`: https://miktex.org/download
- Update `MikTeX`
- Install `SumatraPDF` viewer: https://www.sumatrapdfreader.org/download-free-pdf-viewer
- Path to SumatraPDF: `C:\Program Files\SumatraPDF`
- Install `Intellij IDEA Ultimate`: https://www.jetbrains.com/idea/download/#section=windows
- Activate `Intellij IDEA Ultimate`
- Install `TeXiFy IDEA` plugin: https://plugins.jetbrains.com/plugin/9473-texify-idea
- Clone this repository locally: `https://github.com/kolosovpetro/github-latex-template.git`
- Open `github-latex-template` folder in `Intellij IDEA Ultimate` and configure as follows
    - LaTeX Configuration
      ![LaTeX Configuration](./img/latex_configuration.PNG "LaTeX Configuration")
    - BibTeX Configuration
      ![BibTeX Configuration](./img/bibtex_configuration.PNG "BibTeX Configuration")
- Configure Inverse Search in `Intellij IDEA` for SumatraPDF: `Tools -> LaTeX -> Configure Inverse Search`
- Compile document using `Shift + F10`
