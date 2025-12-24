#!/bin/bash
# MoleculeAI Installation Script

echo "🧬 Installing MoleculeAI..."
echo ""

# Install dependencies
echo "📦 Installing Python packages..."
pip install gradio>=4.0.0
pip install google-generativeai>=0.3.0
pip install python-dotenv>=1.0.0
pip install pandas>=2.0.0
pip install matplotlib>=3.7.0

# Install RDKit with compatible numpy
echo "🧪 Installing RDKit..."
pip install "numpy<2.0"
pip install rdkit-pypi>=2022.9.5

echo ""
echo "✅ Installation complete!"
echo ""
echo "📝 Next steps:"
echo "1. Add your Gemini API key to .env file"
echo "2. Run: python app_gradio.py"
echo ""
