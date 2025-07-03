import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from rdkit_utils import get_rdkit_info
from claude_utils import (
    get_synthesis_info,
    name_to_smiles,
    get_price_estimates,
    get_market_trend
)


import re

def convert_persian_to_english_digits(s):
    persian_to_english = str.maketrans("۰۱۲۳۴۵۶۷۸۹٬", "0123456789,")
    return s.translate(persian_to_english)

import re

def extract_average_price(price_info):
    """
    Extracts the average price from each line of the price_info string.
    Only looks for two English-style numbers per line, computes their average,
    and sums them across all lines.
    """
    total = 0

    for line in price_info.splitlines():
        # Find two numeric values in the line (e.g., 850000 and 1200000)
        match = re.findall(r"(\d[\d,]*)", line)

        if len(match) >= 2:
            try:
                # Remove commas and convert to int
                low = int(match[0].replace(",", ""))
                high = int(match[1].replace(",", ""))
                avg = (low + high) / 2
                total += avg
            except:
                continue  # skip lines with conversion issues

    return round(total)






st.set_page_config(page_title="Molecule Info Explorer", layout="wide")
st.title("🧪 Molecule Info Explorer")

st.markdown("""
Upload a **text file** containing molecule names (one per line) to get their **SMILES**, **molecular properties**, and **industrial synthesis information**.
""")

uploaded_file = st.file_uploader("Upload a .txt file with molecule names", type=["txt"])

# Initialize session state to store results
if "molecule_data" not in st.session_state:
    st.session_state.molecule_data = []

if uploaded_file is not None and not st.session_state.molecule_data:
    molecule_names = [line.decode("utf-8").strip() for line in uploaded_file.readlines()]
    st.markdown(f"📄 Found **{len(molecule_names)}** molecules.")

    results = []

    for name in molecule_names:
        if not name:
            continue
        with st.spinner(f"Processing: {name}"):
            smiles = name_to_smiles(name)
            if not smiles:
                results.append({
                    "Name": name,
                    "SMILES": "Not Found",
                    "Molecular Weight": "N/A",
                    "Formula": "N/A",
                    "Synthetic Route": "N/A",
                    "Raw Materials": "N/A",
                    "Price of Raw Materials": "N/A"
                })
                continue

            chem_info, error = get_rdkit_info(smiles)
            mol_weight = "RDKit error" if error else chem_info.get("Molecular Weight", "N/A")
            formula = "RDKit error" if error else chem_info.get("Molecular Formula", "N/A")

            synthesis = get_synthesis_info(name)
            synthetic_route = synthesis.get("synthetic_route", "N/A")
            raw_materials = ", ".join(synthesis.get("raw_materials", []))
            price_info = get_price_estimates(synthesis.get("raw_materials", []))
            total_price = extract_average_price(price_info)

            results.append({
                "Name": name,
                "SMILES": smiles,
                "Molecular Weight": mol_weight,
                "Formula": formula,
                "Synthetic Route": synthetic_route,
                "Raw Materials": raw_materials,
                "Price of Raw Materials": price_info,
                "Total Price of Raw Materials (Toman)": total_price
            })

    st.session_state.molecule_data = results  # Save to session state

if st.session_state.molecule_data:
    df = pd.DataFrame(st.session_state.molecule_data)
    st.success("✅ Molecule processing complete.")
    st.dataframe(df, use_container_width=True)

    st.markdown("### 📈 Market Analysis Options")

    # Option to analyze a single molecule's market trend
    single_name = st.text_input("🔍 Enter a molecule name to view market trends (2023–2027)")
    if st.button("Show Market Trend for This Molecule"):
        if single_name:
            with st.spinner(f"Fetching market trend for: {single_name}"):
                trend = get_market_trend(single_name)

            if not trend or not trend.get("global_market"):
                st.warning(f"No market data found for {single_name}")
            else:
                years = sorted(trend["global_market"].keys())
                global_vals = [trend["global_market"][y] for y in years]
                me_vals = [trend["middle_east_market"][y] for y in years]

                fig, ax = plt.subplots(figsize=(7, 3), dpi=120)
                ax.plot(years, global_vals, marker='o', linestyle='-', linewidth=2.5, label="Global Market")
                ax.plot(years, me_vals, marker='s', linestyle='--', linewidth=2.5, label="Middle East Market")

                ax.set_title(f"Market Value of {single_name} (2023–2027)", fontsize=14)
                ax.set_ylabel("Market Value (Million USD)")
                ax.set_xlabel("Year")
                ax.grid(True)
                ax.legend()
                st.pyplot(fig)

        else:
            st.warning("Please enter a molecule name.")

else:
    st.info("Please upload a `.txt` file with molecule names.")
