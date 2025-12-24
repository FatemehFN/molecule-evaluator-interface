"""Configuration management for the Molecule Analysis Chatbot."""
import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()

class Config:
    """Application configuration."""

    # API Keys
    GEMINI_API_KEY = os.getenv("GEMINI_API_KEY", "")

    # Paths
    BASE_DIR = Path(__file__).parent
    DATA_DIR = BASE_DIR / "data"
    EQUIPMENT_LIST_PATH = os.getenv(
        "EQUIPMENT_LIST_PATH",
        str(DATA_DIR / "equipment_list.txt")
    )

    # Model Configuration
    GEMINI_MODEL = "gemini-2.0-flash"
    TEMPERATURE = 0.1
    MAX_TOKENS = 8192

    # Rate Limiting
    API_DELAY = 0.1  # seconds between API calls (reduced for faster processing)

    # UI Configuration
    PAGE_TITLE = "MoleculeAI"
    PAGE_ICON = "🧬"
    LAYOUT = "wide"

    @classmethod
    def validate(cls):
        """Validate required configuration."""
        if not cls.GEMINI_API_KEY:
            raise ValueError("GEMINI_API_KEY not found in environment variables")
        return True
