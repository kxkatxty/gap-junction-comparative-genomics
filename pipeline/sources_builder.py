#!/usr/bin/env python3
"""Build the published-literature Sources page for the result sites."""

from __future__ import annotations

from dataclasses import dataclass
from html import escape
from pathlib import Path

from pipeline.common import RESULTS_DIR

OUT = RESULTS_DIR / "sources"


@dataclass(frozen=True)
class Source:
    num: int
    authors: str
    year: int
    title: str
    venue: str
    doi: str | None = None
    url: str | None = None
    mirror_url: str | None = None
    note: str | None = None

    @property
    def primary_url(self) -> str:
        if self.doi:
            return f"https://doi.org/{self.doi}"
        return self.url or "#"

    def citation_html(self) -> str:
        parts = [
            f"{escape(self.authors)} ({self.year}).",
            f"<em>{escape(self.title)}</em>.",
            escape(self.venue) + ".",
        ]
        if self.doi:
            parts.append(
                f'<a href="https://doi.org/{escape(self.doi)}" rel="noopener">https://doi.org/{escape(self.doi)}</a>.'
            )
        elif self.url:
            parts.append(f'<a href="{escape(self.url)}" rel="noopener">{escape(self.url)}</a>.')
        if self.mirror_url:
            parts.append(
                f'Alternate access: <a href="{escape(self.mirror_url)}" rel="noopener">mirror</a>.'
            )
        if self.note:
            parts.append(f"<span class=\"note\">{escape(self.note)}</span>")
        return " ".join(parts)


SECTIONS: tuple[tuple[str, tuple[Source, ...]], ...] = (
    (
        "Reviews & general biology",
        (
            Source(
                1,
                "Bruzzone, R., White, T. W., & Paul, D. L.",
                1996,
                "Connections with connexins: the molecular basis of direct intercellular signaling",
                "European Journal of Biochemistry, 238(1), 1–27",
                doi="10.1111/j.1432-1033.1996.0001q.x",
            ),
            Source(
                2,
                "Dahl, G., & Müller, D. J.",
                2014,
                "Innexin and pannexin channels and their signaling",
                "FEBS Letters, 588(8), 1396–1402",
                doi="10.1016/j.febslet.2014.03.007",
            ),
            Source(
                3,
                "Skerrett, I. M., & Williams, T. C.",
                2017,
                "A structural and functional comparison of gap junction channels composed of connexins and innexins",
                "Developmental Neurobiology, 77(5), 522–547",
                doi="10.1002/dneu.22447",
            ),
            Source(
                9,
                "Scemes, E., Suadicani, S. O., Dahl, G., & Spray, D. C.",
                2007,
                "Connexin and pannexin mediated cell–cell communication",
                "Neuron Glia Biology, 3(3), 199–208",
                doi="10.1017/S1740925X08000069",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2588549/",
            ),
            Source(
                10,
                "Scemes, E., & Spray, D. C.",
                2009,
                "Connexins, pannexins, innexins: novel roles of “hemi-channels”",
                "Pflügers Archiv, 457(6), 1207–1216",
                doi="10.1007/s00424-008-0591-5",
            ),
            Source(
                11,
                "Wei, C. J., Xu, X., & Lo, C. W.",
                2004,
                "Connexins and cell signaling in development and disease",
                "Annual Review of Cell and Developmental Biology, 20, 811–838",
                doi="10.1146/annurev.cellbio.19.111301.144309",
            ),
            Source(
                12,
                "Goodenough, D. A., Goliger, J. A., & Paul, D. L.",
                1996,
                "Connexins, connexons, and intercellular communication",
                "Annual Review of Biochemistry, 65, 475–502",
                doi="10.1146/annurev.bi.65.070196.002355",
            ),
            Source(
                13,
                "Sáez, J. C., Berthoud, V. M., Branes, M. C., Martinez, A. D., & Beyer, E. C.",
                2003,
                "Plasma membrane channels formed by connexins: their regulation and functions",
                "Physiological Reviews, 83(4), 1359–1400",
                doi="10.1152/physrev.00007.2003",
            ),
            Source(
                14,
                "Maeda, S., & Tsukihara, T.",
                2011,
                "Structure of the gap junction channel and its implications for its biological functions",
                "Cellular and Molecular Life Sciences, 68(7), 1115–1129",
                doi="10.1007/s00018-010-0551-z",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3071404/",
            ),
            Source(
                15,
                "Falk, M. M., Kells, R. M., Berthoud, V. M., Kanaporis, G., & Beyer, E. C.",
                2013,
                "Proteins and mechanisms regulating gap-junction assembly, internalization, and degradation",
                "Physiology, 28(2), 93–107",
                doi="10.1152/physiol.00038.2012",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC3642965/",
            ),
            Source(
                16,
                "Esseltine, J. L., & Laird, D. W.",
                2017,
                "Gap junction gene and protein families: connexins, innexins, and pannexins",
                "Biochimica et Biophysica Acta — Biomembranes, 1860(1), 5–8",
                doi="10.1016/j.bbamem.2017.05.016",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC5704981/",
            ),
            Source(
                21,
                "Peng, B., Xu, C., Wang, S., Zhang, Y., & Li, W.",
                2022,
                "The role of connexin hemichannels in inflammatory diseases",
                "Biology, 11(2), 237",
                doi="10.3390/biology11020237",
            ),
            Source(
                20,
                "Hanner, N., et al.",
                2010,
                "Connexins and the kidney",
                "American Journal of Physiology — Renal Physiology, 298(3), F572–F582",
                doi="10.1152/ajpregu.00808.2009",
            ),
        ),
    ),
    (
        "Evolution & comparative context",
        (
            Source(
                4,
                "Welzel, L. F., & Schuster, S. C.",
                2022,
                "Connexins evolved after early chordates lost innexin diversity",
                "eLife, 11, e74422",
                doi="10.7554/eLife.74422",
            ),
            Source(
                29,
                "Sohl, G., & Willecke, K.",
                2003,
                "Gap junctions and the connexin gene family",
                "Cell Communication & Adhesion, 10(4–6), 173–180",
                doi="10.1080/cac.10.4-6.173.180",
            ),
            Source(
                30,
                "Starich, T. A., Miller, A., Nguyen, R. L., Davis, M. W., & Hall, D. H.",
                2001,
                "The innexin family of invertebrate gap junction proteins",
                "Cell Communication & Adhesion, 8(4–6), 303–308",
                doi="10.3109/15419060109080744",
            ),
        ),
    ),
    (
        "High-resolution structures",
        (
            Source(
                5,
                "Oshima, C., Tani, K., Hiroaki, Y., Fujiyoshi, Y., & Doi, T.",
                2019,
                "Cryo-EM structures of undocked innexin-6 hemichannels in phospholipids",
                "Science Advances, 5(8), eaax3157",
                doi="10.1126/sciadv.aax3157",
            ),
            Source(
                6,
                "Deng, W., Lu, W., Du, J., Ruan, Z., et al.",
                2020,
                "Cryo-EM structure of human heptameric Pannexin 1 channel",
                "Cell Research, 30, 446–448",
                doi="10.1038/s41422-020-0298-5",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC7185789/",
            ),
            Source(
                7,
                "Michalski, K., Syrjanen, J. L., Henze, E., Kumpf, J., Furukawa, H., & Kawate, T.",
                2020,
                "The Cryo-EM structure of pannexin 1 reveals unique motifs for ion selection and inhibition",
                "eLife, 9, e54670",
                doi="10.7554/eLife.54670",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC7108861/",
            ),
            Source(
                8,
                "Bennett, B. C., Harris, A. L., et al.",
                2024,
                "Connexin gap junction channels and hemichannels: insights from high-resolution structures",
                "Biology, 13(5), 298",
                doi="10.3390/biology13050298",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC11120952/",
            ),
        ),
    ),
    (
        "Innexin biology (invertebrate models)",
        (
            Source(
                23,
                "Güiza, J., Barría, I., Sáez, J. C., & Vega, J. L.",
                2018,
                "Innexins: expression, regulation, and functions",
                "Frontiers in Physiology, 9, 1414",
                doi="10.3389/fphys.2018.01414",
                url="https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2018.01414/full",
            ),
            Source(
                24,
                "Güiza, J., Barría, I., Sáez, J. C., & Vega, J. L.",
                2018,
                "Innexins: expression, regulation, and functions",
                "Frontiers in Physiology, 9, 1414 (PMC mirror)",
                url="https://pmc.ncbi.nlm.nih.gov/articles/PMC6193117/",
                note="Same review as #23; alternate index via PubMed Central.",
            ),
            Source(
                25,
                "Bauer, R., Lehmann, C., Martini, J., & Eckardt, F.",
                2004,
                "Gap junction channel protein innexin 2 is essential for epithelial morphogenesis in the Drosophila embryo",
                "Molecular Biology of the Cell, 15(6), 2992–3004",
                doi="10.1091/mbc.E04-01-0056",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC420120/",
            ),
            Source(
                26,
                "Giuliani, F., Giuliani, G., Bauer, R., & Rabouille, C.",
                2013,
                "Innexin 3, a new gene required for dorsal closure in Drosophila embryo",
                "PLOS ONE, 8(7), e69212",
                doi="10.1371/journal.pone.0069212",
            ),
            Source(
                27,
                "Lehmann, C., Leitner, J., Altenhein, B., & Bauer, R.",
                2006,
                "Heteromerization of innexin gap junction proteins regulates epithelial tissue organization in Drosophila",
                "Molecular Biology of the Cell, 17(4), 1676–1685",
                doi="10.1091/mbc.E05-11-1059",
                mirror_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC1415333/",
            ),
            Source(
                28,
                "Kelley, T. W., et al.",
                2022,
                "Innexin function dictates the spatial relationship between distal somatic cells in the Caenorhabditis elegans gonad without impacting the germline stem cell pool",
                "eLife, 11, e74955",
                doi="10.7554/eLife.74955",
            ),
            Source(
                22,
                "Walker, D. S., et al.",
                2020,
                "Distinct roles for innexin gap junctions and hemichannels in mechanosensation",
                "eLife, 9, e50597",
                doi="10.7554/eLife.50597",
            ),
        ),
    ),
    (
        "Books & reference works",
        (
            Source(
                18,
                "Harris, A. L. (Ed.)",
                2007,
                "Connexins: a guide",
                "Methods in Molecular Biology, vol. 407. Humana Press",
                doi="10.1007/978-1-59745-489-6",
                url="https://link.springer.com/book/10.1007/978-1-59745-489-6",
                note="Methods in Molecular Biology volume; also indexed on Google Books.",
            ),
            Source(
                19,
                "Dhein, S., Mohr, F. W., & Delmar, M. (Eds.)",
                2006,
                "Cardiovascular gap junctions",
                "Advances in Cardiology, vol. 42. Karger",
                url="https://www.karger.com/Book/Home/255842",
                note="ISBN 978-3-8055-8077-0.",
            ),
        ),
    ),
    (
        "Web & archived resources",
        (
            Source(
                17,
                "Wikipedia contributors",
                2024,
                "Gap junction",
                "Wikipedia, The Free Encyclopedia",
                url="https://en.wikipedia.org/wiki/Gap_junction",
            ),
            Source(
                31,
                "Stagg, M. A., & Fletcher, W. H.",
                1990,
                "The relationship between connexins and gap junction intercellular communication",
                "Endocrine Reviews, 11(2), 302–325",
                doi="10.1210/edrv-11-2-302",
                mirror_url="https://web.archive.org/web/20180727114153id_/https://academic.oup.com/edrv/article-pdf/11/2/302/4912342/11-2-302.pdf",
                note="Original Silverchair/OUP PDF link was truncated in the source list; this is the best-matching archived copy.",
            ),
        ),
    ),
)


def _count_sources() -> int:
    seen: set[int] = set()
    for _, sources in SECTIONS:
        for src in sources:
            seen.add(src.num)
    return len(seen)


def build(out_dir: Path = OUT) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    n_sources = _count_sources()

    sections_html: list[str] = []
    for title, sources in SECTIONS:
        items = "\n".join(
            f'      <li id="source-{src.num}" value="{src.num}"><p>{src.citation_html()}</p></li>'
            for src in sources
        )
        sections_html.append(
            f"""    <section class="group">
      <h2>{escape(title)}</h2>
      <ol class="refs">
{items}
      </ol>
    </section>"""
        )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Sources — gap-junction comparative genomics</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600;700&family=Sora:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root {{
  --ink:#12202a; --muted:#4d6570; --sea:#0f766e; --line:rgba(18,32,42,.12);
  --card:#fff; --foam:#f4f7f5;
}}
* {{ box-sizing:border-box }}
body {{
  margin:0; color:var(--ink); font-family:Sora,sans-serif; line-height:1.65;
  background:linear-gradient(180deg,#e8f1ee 0%, var(--foam) 45%, #eef3f1 100%);
}}
a {{ color:var(--sea) }}
.wrap {{ width:min(820px, calc(100% - 2rem)); margin:0 auto; padding:2.2rem 0 4rem }}
.hero {{
  padding:2rem 0 1.2rem; border-bottom:1px solid var(--line); margin-bottom:1.6rem;
}}
.hero h1 {{
  font-family:Fraunces,serif; font-size:clamp(1.8rem,4vw,2.4rem);
  margin:0 0 .45rem; letter-spacing:-.02em;
}}
.hero p {{ margin:0; color:var(--muted); max-width:58ch }}
.nav {{ margin:1rem 0 0; font-size:.88rem }}
.nav a {{ margin-right:.85rem; font-weight:600; text-decoration:none }}
.group {{ margin:1.8rem 0 }}
.group h2 {{
  font-family:Fraunces,serif; font-size:1.15rem; margin:0 0 .75rem; color:var(--ink);
}}
.refs {{
  list-style:none; margin:0; padding:0; display:grid; gap:.85rem;
}}
.refs li {{
  background:var(--card); border:1px solid var(--line); border-radius:14px;
  padding:.85rem 1rem; box-shadow:0 6px 18px rgba(18,32,42,.04);
}}
.refs li::before {{
  content: attr(value);
  display:inline-block; min-width:1.6rem; font-weight:700; color:var(--sea); margin-right:.35rem;
}}
.refs p {{ margin:0; font-size:.92rem }}
.note {{ display:block; margin-top:.35rem; font-size:.82rem; color:var(--muted) }}
.footer {{ margin-top:2rem; color:var(--muted); font-size:.88rem }}
.footer a {{ margin-right:.85rem }}
</style>
</head>
<body>
<main class="wrap">
  <header class="hero">
    <h1>Sources</h1>
    <p>
      Published literature and web references used for background on connexins, innexins, and pannexins
      on this site. {n_sources} numbered entries from the thesis reading list.
    </p>
    <p class="nav">
      <a href="../gap_junction_path/index.html">← Main path</a>
      <a href="../inx_vs_cnx/index.html">Innexin × connexin</a>
      <a href="../phylogenetic_story/index.html">Phylogeny story</a>
    </p>
  </header>

{chr(10).join(sections_html)}

  <footer class="footer">
    <p>Figures and sequence comparisons on other pages come from this repository’s pipeline outputs; this page lists external literature only.</p>
    <p>
      <a href="../gap_junction_path/index.html">Gap-junction path</a>
      <a href="../../../pannexin/project/results/pannexin_path/index.html">Pannexin path</a>
      <a href="../../../pannexin/project/results/panx_vs_inx/index.html">Pannexin × innexin</a>
    </p>
  </footer>
</main>
</body>
</html>
"""
    (out_dir / "index.html").write_text(html, encoding="utf-8")
    print(f"Wrote {out_dir / 'index.html'} ({n_sources} sources)")


if __name__ == "__main__":
    build()
