import anthropic
import time
import json
import re

# Claude API key – stored here for development use (move to env var for production)
CLAUDE_API_KEY = "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"

client = anthropic.Anthropic(api_key=CLAUDE_API_KEY)

def get_synthesis_info(molecule_name):
    """
    Uses Claude to get industrial synthesis route and raw materials for a molecule.
    Returns a dictionary with keys: 'synthetic_route' and 'raw_materials'.
    """
    try:
        time.sleep(2.5)  # Slight pause to avoid potential overuse/rate-limiting

        response = client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=512,
            temperature=0.1,
            system="You are a chemistry assistant. Always respond with accurate and validated industrial synthesis steps and raw materials in strict JSON format.",
            messages=[
                {
                    "role": "user",
                    "content": f"""Provide the correct industrial synthesis route and raw materials for {molecule_name}. Use this JSON format exactly:

{{
  "synthetic_route": [
    "1. [Reaction Step 1]",
    "2. [Reaction Step 2]",
    "...",
    "n. [Final Step]"
  ],
  "raw_materials": ["Material 1", "Material 2", "..."]
}}

Only return the JSON. Do not include any markdown or comments."""
                }
            ]
        )

        # Attempt to parse JSON-like string output
        import json
        import re

        raw = response.content[0].text.strip()
        json_text = re.search(r"\{.*\}", raw, re.DOTALL)
        if json_text:
            return json.loads(json_text.group())
        else:
            return {
                "synthetic_route": ["Could not parse response."],
                "raw_materials": [],
                "error": "Invalid JSON format returned by Claude"
            }

    except Exception as e:
        return {
            "synthetic_route": ["Not found (error)"],
            "raw_materials": ["Not available"],
            "error": str(e)
        }



def name_to_smiles(molecule_name):
    try:
        response = client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=100,
            temperature=0.1,
            system="You are a chemistry assistant that only returns SMILES strings.",
            messages=[
                {"role": "user", "content": f"Give only the SMILES string for this molecule: {molecule_name}"}
            ]
        )
        text = response.content[0].text.strip()
        # optionally sanitize or validate
        return text
    except Exception as e:
        return None
    



def get_price_estimates(raw_materials):
    try:
        materials_list = "\n".join(f"- {m}" for m in raw_materials)
        response = client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=1024,
            temperature=0.2,
            system="You are a chemical market analyst fluent in Persian (Farsi). You search Iranian supplier websites for average price estimates in toman or rials.",
            messages=[
                {
                    "role": "user",
                    "content": f"""برای هر یک از مواد اولیه زیر، لطفاً میانگین بازه قیمت آن را از وب‌سایت‌های فروشنده‌های ایرانی مانند شیمی دشت، کیان تجهیز شیمی، بهان‌سر، مرک ایران و دیگر فروشگاه‌های آنلاین ایرانی پیدا کن.

مواد اولیه:
{materials_list}

فرمت پاسخ:
نام ماده: بازه قیمت تقریبی (نام شرکت ۱، نام شرکت ۲، ...)

مثال:
اسید استیک: ۵۰۰٬۰۰۰–۷۰۰٬۰۰۰ تومان در هر کیلو (شیمی دشت، کیان تجهیز)

فقط فهرست قیمت را بده. از نوشتن توضیح اضافی خودداری کن."""
                }
            ]
        )

        return response.content[0].text.strip()
    except Exception as e:
        return f"⚠️ قیمت مواد در دسترس نیست (خطا: {e})"





def get_market_trend(material):
    try:
        response = client.messages.create(
            model="claude-3-5-sonnet-20241022",
            max_tokens=1024,
            temperature=0.1,
            system="You are an expert chemical market analyst. Estimate the market value (in million USD) of a chemical material in the global and Middle East markets for 2023 to 2027. Return a clean JSON.",
            messages=[
                {
                    "role": "user",
                    "content": f"""Give the estimated market value for {material} from 2023 to 2027 in million USD.
Break it down into:
- global_market: {{year: value, ...}}
- middle_east_market: {{year: value, ...}}
Respond in this format:
{{
  \"global_market\": {{\"2023\": 1.2, \"2024\": 1.4, ...}},
  \"middle_east_market\": {{\"2023\": 0.3, \"2024\": 0.4, ...}}
}}"""
                }
            ]
        )

        raw = response.content[0].text.strip()
        json_text = re.search(r"\{.*\}", raw, re.DOTALL)
        if json_text:
            return json.loads(json_text.group())
        else:
            return {}

    except Exception as e:
        return {}