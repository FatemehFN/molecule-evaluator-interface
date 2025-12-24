"""Base agent class for the molecule analysis system."""
import os
import time
import json
import re
from typing import Dict, Any, Optional
from abc import ABC, abstractmethod
import requests
from dotenv import load_dotenv

# Load environment variables
from pathlib import Path
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(dotenv_path=env_path, override=True)

from config import Config


class BaseAgent(ABC):
    """Abstract base class for all agents using REST API."""

    def __init__(self, name: str, role: str):
        """Initialize the agent.

        Args:
            name: Agent name
            role: Agent role/specialty
        """
        self.name = name
        self.role = role
        self.api_key = Config.GEMINI_API_KEY
        self.model_name = Config.GEMINI_MODEL
        self.api_url = f"https://generativelanguage.googleapis.com/v1beta/models/{self.model_name}:generateContent"
        self.timeout = 60  # 60 second timeout for API calls

    def _call_api(self, prompt: str, system_instruction: Optional[str] = None, timeout: Optional[int] = None) -> str:
        """Make API call to Gemini using REST API.

        Args:
            prompt: User prompt
            system_instruction: System instruction for the model
            timeout: Timeout in seconds (uses default if not specified)

        Returns:
            Generated text response
        """
        if timeout is None:
            timeout = self.timeout
            
        time.sleep(Config.API_DELAY)

        # Build the request payload
        full_prompt = prompt
        if system_instruction:
            full_prompt = f"{system_instruction}\n\n{prompt}"

        payload = {
            "contents": [{
                "parts": [{"text": full_prompt}]
            }],
            "generationConfig": {
                "temperature": Config.TEMPERATURE,
                "maxOutputTokens": Config.MAX_TOKENS
            }
        }

        try:
            response = requests.post(
                f"{self.api_url}?key={self.api_key}",
                json=payload,
                timeout=timeout,
                headers={'Content-Type': 'application/json'}
            )

            if response.status_code == 503:
                return "Error: Gemini API overloaded (503). Please try again in a moment."
            
            if response.status_code != 200:
                return f"Error: API request failed with status {response.status_code}: {response.text[:200]}"

            result = response.json()
            
            # Extract text from response
            if 'candidates' in result and result['candidates']:
                candidate = result['candidates'][0]
                if 'content' in candidate and candidate['content'].get('parts'):
                    return candidate['content']['parts'][0]['text'].strip()
                else:
                    finish_reason = candidate.get('finishReason', 'UNKNOWN')
                    return f"Error: Content blocked. Reason: {finish_reason}"
            else:
                return f"Error: Unexpected API response structure: {result}"

        except requests.Timeout:
            return f"Error: API call timed out after {timeout} seconds"
        except requests.ConnectionError as e:
            return f"Error: Connection failed: {str(e)}"
        except Exception as e:
            return f"Error: {str(e)}"

    def _extract_json(self, text: str) -> Optional[Dict[str, Any]]:
        """Extract JSON from text response.

        Args:
            text: Text containing JSON

        Returns:
            Parsed JSON dict or None
        """
        json_match = re.search(r'\{.*\}', text, re.DOTALL)
        if json_match:
            try:
                return json.loads(json_match.group())
            except json.JSONDecodeError:
                return None
        return None

    @abstractmethod
    def process(self, *args, **kwargs) -> Any:
        """Process the agent's task. Must be implemented by subclasses."""
        pass
