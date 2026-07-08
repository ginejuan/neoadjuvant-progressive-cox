#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
Supplementary code S2 — Source-document data-recovery extractors
Manuscript: "Progressive Cox models expose confounding by indication in
neoadjuvant chemotherapy ..." (Clinical and Translational Oncology)
================================================================================

Covariate values missing in the clinical database were recovered from the
original source documents (Supplementary Figure S1). This file documents the
text-extraction functions applied to the plain text of the source PDFs
(pathology reports and oncology anamnesis/evolution notes), obtained with
`pdftotext` (Poppler).

NOTE ON REPRODUCIBILITY AND PRIVACY
-----------------------------------
These functions operate on the ORIGINAL clinical PDFs, which contain personal
health information and therefore CANNOT be shared. They are provided for
methodological transparency only. Each automated value was validated against
the patients with the value already recorded, with the concordance reported in
the manuscript (grade 92% exact; Ki67 r=0.93 and 98% agreement at the 20%
threshold; age mean absolute error 0.06 years; radiotherapy manually confirmed
against the radiotherapy register).

Dependencies: pdftotext (Poppler) for text extraction; Python standard library.
================================================================================
"""
import re
import subprocess


def pdf_text(path, max_pages=8):
    """Plain text of the first `max_pages` pages of a PDF (Poppler pdftotext)."""
    try:
        out = subprocess.run(["pdftotext", "-l", str(max_pages), path, "-"],
                             capture_output=True, timeout=15)
        return out.stdout.decode("utf-8", "ignore")
    except Exception:
        return ""


# --- Histological grade (Nottingham / Bloom-Richardson) ---------------------
def extract_grade(text):
    u = text.upper()
    m = re.search(r"GRADO\s+(?:HISTOL[OÓ]GICO\s+)?(?:DE\s+)?"
                  r"(?:NOTTINGHAM\s+|SBR\s+|ELSTON[- ]ELLIS\s+)?(III|II|I)\b", u)
    if m:
        return {"I": 1, "II": 2, "III": 3}[m.group(1)]
    m = re.search(r"GRADO\s+(?:HISTOL[OÓ]GICO\s+)?(?:DE\s+)?(?:NOTTINGHAM\s+|SBR\s+)?([123])\b", u)
    if m:
        return int(m.group(1))
    m = re.search(r"\bG\s*([123])\b", u)
    if m:
        return int(m.group(1))
    if "POBREMENTE DIFERENCIADO" in u or "POCO DIFERENCIADO" in u:
        return 3
    if "MODERADAMENTE DIFERENCIADO" in u:
        return 2
    if "BIEN DIFERENCIADO" in u:
        return 1
    return None


# --- Ki67 proliferation index (%) -------------------------------------------
def extract_ki67(text):
    u = text.upper()
    patterns = [
        (r"KI\s*-?\s*67[^\n%]{0,40}?(\d{1,3})\s*%?\s*-\s*(\d{1,3})\s*%", "range"),
        (r"KI\s*-?\s*67:?\s*POSITIVO\s+EN\s+(\d{1,3})\s*%", "one"),
        (r"KI\s*-?\s*67[^\n]{0,30}?(?:DEL|EN|:|=)\s*(\d{1,3})\s*%", "one"),
        (r"KI\s*-?\s*67[:\s]{0,6}(\d{1,3})\s*%", "one"),
        (r"(?:[ÍI]NDICE\s+DE\s+PROLIFERACI[ÓO]N|PROLIFERACI[ÓO]N)[^\n]{0,30}?"
         r"KI\s*-?\s*67[^\n]{0,20}?(\d{1,3})\s*%", "one"),
        (r"MIB\s*-?\s*1[^\n]{0,20}?(\d{1,3})\s*%", "one"),
        (r"KI\s*-?\s*67\s*\)?\s*,\s*(\d{1,3})\s*%", "one"),      # "(KI67), 12%"
        (r"KI\s*-?\s*67\s*\)\s*:?\s*(\d{1,3})\s*%", "one"),       # "(KI67): 12%"
    ]
    for pat, kind in patterns:
        m = re.search(pat, u)
        if m:
            if kind == "range":
                a, b = int(m.group(1)), int(m.group(2))
                if a <= 100 and b <= 100:
                    return round((a + b) / 2)
            else:
                v = int(m.group(1))
                if 0 <= v <= 100:
                    return v
    return None


# --- Date of birth -> age at diagnosis --------------------------------------
def extract_dob(text):
    """Date of birth (dd/mm/yyyy) from the oncology anamnesis header."""
    m = re.search(r"F\.?\s*NACIMIENTO:?\s*(\d{2})/(\d{2})/(\d{4})", text.upper())
    if m:
        return (int(m.group(1)), int(m.group(2)), int(m.group(3)))  # (day, month, year)
    return None
# age at diagnosis = (biopsy_date - date_of_birth).days / 365.25


# --- Adjuvant radiotherapy (free-text classifier; requires manual review) ---
RT_RECEIVED = [r"TRAS\s+(?:LA\s+)?RADIOTERAPIA", r"FINALIZ[AÓ][A-Z]*\s+(?:LA\s+)?RADIOTERAPIA",
               r"HA\s+FINALIZADO\s+LA\s+RADIOTERAPIA", r"FIN\s+DE\s+(?:LA\s+)?RADIOTERAPIA",
               r"LLEVA[^.]{0,25}CON\s+RADIOTERAPIA", r"RADIODERMITIS",
               r"RECIBI[ÓO][^.]{0,25}RADIOTERAPIA", r"\d+\s*GY", r"SESIONES\s+DE\s+RADIOTERAPIA"]
RT_NEG = [r"NO\s+(?:ES\s+|SE\s+|PRECISA|CANDIDATA|TRIBUTARIA)[^.]{0,25}RADIOTERAPIA",
          r"RECHAZA[^.]{0,20}RADIOTERAPIA", r"CONTRAINDICAD[AO][^.]{0,20}RADIOTERAPIA"]


def classify_radiotherapy(text):
    """Returns 'received', 'no' or 'uncertain'. Free-text; every case was
    manually confirmed against the radiotherapy register before integration."""
    u = re.sub(r"\s+", " ", text.upper())
    if any(re.search(p, u) for p in RT_NEG):
        return "no"
    if any(re.search(p, u) for p in RT_RECEIVED):
        return "received"
    return "uncertain"
