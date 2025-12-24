#!/bin/bash
# MoleculeAI - FAST Mode Launcher

echo "⚡ Starting MoleculeAI FAST Mode..."
echo ""

# Kill any existing processes
echo "🔄 Checking for existing processes..."
lsof -ti:7860 | xargs kill -9 2>/dev/null
sleep 1

# Start the fast app
echo "🚀 Launching FAST mode..."
echo "📍 Access at: http://localhost:7860"
echo ""
echo "⚡ Speed: 3-4x faster than standard mode"
echo "📊 Only 2 API calls per molecule"
echo ""

python app_fast.py
