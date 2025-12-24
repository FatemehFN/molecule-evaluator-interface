#!/bin/bash
# MoleculeAI - REST API Launcher

echo "🚀 Starting MoleculeAI REST API..."
echo ""

# Kill any existing processes
lsof -ti:8000 | xargs kill -9 2>/dev/null
sleep 1

echo "📍 API will run at: http://localhost:8000"
echo "📖 API docs at: http://localhost:8000/docs"
echo "🌐 Web interface: open index.html in browser"
echo ""
echo "⚡ Speed: 3-4x faster (2 API calls per molecule)"
echo "📊 Features: Single analysis, batch processing, AI chat"
echo ""

python api.py
