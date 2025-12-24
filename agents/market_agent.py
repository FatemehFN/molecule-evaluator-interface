"""Market Agent for pricing and market analysis."""
from typing import Dict, List
from .base_agent import BaseAgent


class MarketAgent(BaseAgent):
    """Agent specialized in market analysis and pricing."""

    def __init__(self):
        super().__init__(
            name="MarketAgent",
            role="Chemical market analysis and pricing specialist"
        )

    def get_raw_material_prices(self, raw_materials: List[str]) -> str:
        """Get price estimates for raw materials from Iranian suppliers.

        Args:
            raw_materials: List of raw material names

        Returns:
            Formatted price information in Persian
        """
        if not raw_materials:
            return "No raw materials specified"

        materials_list = "\n".join(f"- {m}" for m in raw_materials)

        system_instruction = (
            "You are a chemical market analyst fluent in Persian (Farsi). "
            "You provide realistic price estimates for chemical raw materials "
            "based on Iranian market data."
        )

        prompt = f"""برای هر یک از مواد اولیه زیر، میانگین قیمت بازار ایران را تخمین بزن:

{materials_list}

فرمت پاسخ:
• نام ماده: [بازه قیمت] تومان به ازای هر کیلوگرم

مثال:
• اسید استیک: ۵۰۰,۰۰۰–۷۰۰,۰۰۰ تومان

فقط لیست قیمت‌ها را بنویس. بدون توضیحات اضافی."""

        return self._call_api(prompt, system_instruction)

    def get_market_trend(self, molecule_name: str) -> Dict:
        """Get market trend data for a molecule.

        Args:
            molecule_name: Name of the molecule

        Returns:
            Dictionary with market data
        """
        system_instruction = (
            "You are a chemical market analyst. Provide realistic market "
            "value estimates in million USD. Return only valid JSON."
        )

        prompt = f"""Estimate the market value for {molecule_name} from 2023 to 2027.

Return ONLY valid JSON in this format:
{{
  "global_market": {{"2023": value, "2024": value, "2025": value, "2026": value, "2027": value}},
  "middle_east_market": {{"2023": value, "2024": value, "2025": value, "2026": value, "2027": value}}
}}

Values in million USD. Be realistic based on the molecule's industrial importance."""

        response = self._call_api(prompt, system_instruction)
        data = self._extract_json(response)

        return data if data else {}

    def calculate_total_price(self, price_text: str) -> int:
        """Calculate total average price from Persian price text.

        Args:
            price_text: Price information in Persian

        Returns:
            Total average price in Toman
        """
        import re

        # Convert Persian digits to English
        persian_to_english = str.maketrans("۰۱۲۳۴۵۶۷۸۹٬", "0123456789,")
        normalized = price_text.translate(persian_to_english)

        total = 0
        # Find all number pairs (ranges)
        matches = re.findall(r"(\d[\d,]*)\s*[-–—]\s*(\d[\d,]*)", normalized)

        for low_str, high_str in matches:
            try:
                low = int(low_str.replace(",", ""))
                high = int(high_str.replace(",", ""))
                total += (low + high) / 2
            except ValueError:
                continue

        return round(total)

    def process(self, raw_materials: List[str], molecule_name: str = None) -> Dict:
        """Process complete market analysis.

        Args:
            raw_materials: List of raw materials
            molecule_name: Optional molecule name for trend analysis

        Returns:
            Dictionary with market analysis
        """
        price_info = self.get_raw_material_prices(raw_materials)
        total_price = self.calculate_total_price(price_info)

        result = {
            "price_info": price_info,
            "total_price_toman": total_price
        }

        if molecule_name:
            result["market_trend"] = self.get_market_trend(molecule_name)

        return result
