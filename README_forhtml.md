# Ullmann C–N Coupling Chemical Space Explorer

Static GitHub Pages dashboard generated from the supplied amine library,
aryl-bromide library, and `all_products_for_ullmann1_2.xlsx`.

## What it does
- accepts amine and aryl bromide SMILES;
- uses RDKit.js in the browser to parse/canonicalize the SMILES;
- calculates Morgan fingerprints (radius 2, 2048 bits);
- uses exact canonical-SMILES matching when available;
- otherwise finds the closest library molecule by Tanimoto similarity;
- shows the user structures and match information;
- places the matched substrate pair on the two-descriptor threshold plot;
- projects it into the PCA chemical space using the fitted scaler/PCA parameters.

RDKit.js is the official JavaScript/WebAssembly distribution of RDKit.

## Test locally
From this folder:

```bash
python -m http.server 8000
```

Then open `http://localhost:8000`.

Do not open `index.html` by double-clicking because browser security can block
loading the JSON data file.

## GitHub Pages
Push the contents of this folder to the desired GitHub repository, then enable
GitHub Pages in Repository Settings → Pages.

## Current data note
The PCA reference plot in this build uses the rows available in
`all_products_for_ullmann1_2.xlsx`. If the notebook's `data_df` contains a
larger enumerated chemical-space background, provide that source file and it
can be added without changing the user-input workflow.
