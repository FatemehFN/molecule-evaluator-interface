"""MoleculeAI - FastAPI REST API.

Fast, professional REST API for molecule analysis.
Perfect for internal use, integrations, and building custom frontends.
"""

from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Optional
import uvicorn
import google.generativeai as genai
from datetime import datetime
import uuid

from config import Config
from agents import UnifiedAgent, ChatAgent


# Initialize configuration
try:
    Config.validate()
    genai.configure(api_key=Config.GEMINI_API_KEY)
except ValueError as e:
    print(f"⚠️ Configuration Error: {e}")
    exit(1)


# Initialize FastAPI
app = FastAPI(
    title="MoleculeAI API",
    description="Fast REST API for comprehensive molecule analysis",
    version="2.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# Enable CORS for web access
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# Initialize agents
unified_agent = UnifiedAgent()
chat_agent = ChatAgent()

# Store analysis results (in-memory for now)
analysis_cache = {}


# Request/Response Models
class MoleculeRequest(BaseModel):
    name: str


class BatchRequest(BaseModel):
    molecules: List[str]


class ChatRequest(BaseModel):
    message: str
    session_id: Optional[str] = None


class AnalysisResponse(BaseModel):
    molecule_name: str
    status: str
    smiles: Optional[str] = None
    properties: Optional[Dict] = None
    synthesis: Optional[Dict] = None
    market: Optional[Dict] = None
    feasibility: Optional[Dict] = None
    applications: Optional[str] = None
    analysis_time: Optional[float] = None


class BatchResponse(BaseModel):
    total: int
    successful: int
    failed: int
    results: List[AnalysisResponse]
    total_time: float


# Root endpoint
@app.get("/")
async def root():
    """API information."""
    return {
        "name": "MoleculeAI API",
        "version": "2.0.0",
        "status": "running",
        "endpoints": {
            "analyze": "/api/analyze",
            "batch": "/api/batch",
            "chat": "/api/chat",
            "docs": "/docs"
        },
        "features": [
            "Fast molecule analysis (2 API calls per molecule)",
            "Batch processing",
            "AI chat assistant",
            "RESTful API"
        ]
    }


@app.get("/health")
async def health_check():
    """Health check endpoint."""
    return {"status": "healthy", "timestamp": datetime.now().isoformat()}


@app.post("/api/analyze", response_model=AnalysisResponse)
async def analyze_molecule(request: MoleculeRequest):
    """Analyze a single molecule.

    Fast analysis using unified agent (only 2 API calls).

    Example:
        POST /api/analyze
        {"name": "Ascorbic acid"}
    """
    start_time = datetime.now()

    try:
        # Process molecule
        result = unified_agent.process(request.name)

        # Calculate processing time
        end_time = datetime.now()
        analysis_time = (end_time - start_time).total_seconds()

        # Store in cache
        cache_id = str(uuid.uuid4())
        analysis_cache[cache_id] = result

        return AnalysisResponse(
            molecule_name=result.get("name"),
            status=result.get("status"),
            smiles=result.get("smiles"),
            properties=result.get("properties"),
            synthesis=result.get("synthesis"),
            market=result.get("market"),
            feasibility=result.get("feasibility"),
            applications=result.get("applications"),
            analysis_time=analysis_time
        )

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/api/batch", response_model=BatchResponse)
async def analyze_batch(request: BatchRequest):
    """Analyze multiple molecules in batch.

    Example:
        POST /api/batch
        {"molecules": ["Ascorbic acid", "Retinol", "Niacinamide"]}
    """
    start_time = datetime.now()
    results = []
    successful = 0
    failed = 0

    for molecule_name in request.molecules:
        try:
            result = unified_agent.process(molecule_name)

            if result.get("status") == "success":
                successful += 1
            else:
                failed += 1

            results.append(AnalysisResponse(
                molecule_name=result.get("name"),
                status=result.get("status"),
                smiles=result.get("smiles"),
                properties=result.get("properties"),
                synthesis=result.get("synthesis"),
                market=result.get("market"),
                feasibility=result.get("feasibility"),
                applications=result.get("applications")
            ))

        except Exception as e:
            failed += 1
            results.append(AnalysisResponse(
                molecule_name=molecule_name,
                status="failed",
                applications=f"Error: {str(e)}"
            ))

    end_time = datetime.now()
    total_time = (end_time - start_time).total_seconds()

    return BatchResponse(
        total=len(request.molecules),
        successful=successful,
        failed=failed,
        results=results,
        total_time=total_time
    )


@app.post("/api/chat")
async def chat(request: ChatRequest):
    """Chat with AI about analyzed molecules.

    Example:
        POST /api/chat
        {"message": "Which molecule is cheapest?", "session_id": "optional-id"}
    """
    try:
        # Get molecules from cache
        molecules = list(analysis_cache.values())

        if not molecules:
            return {
                "response": "Please analyze some molecules first before asking questions.",
                "session_id": request.session_id or str(uuid.uuid4())
            }

        # Simple chat (single message for now)
        messages = [{"role": "user", "content": request.message}]
        response = chat_agent.process(messages, molecules)

        return {
            "response": response,
            "session_id": request.session_id or str(uuid.uuid4()),
            "molecules_in_context": len(molecules)
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/api/molecules")
async def list_analyzed_molecules():
    """Get list of all analyzed molecules."""
    molecules = [
        {
            "name": m.get("name"),
            "status": m.get("status"),
            "smiles": m.get("smiles")
        }
        for m in analysis_cache.values()
    ]

    return {
        "total": len(molecules),
        "molecules": molecules
    }


@app.delete("/api/cache")
async def clear_cache():
    """Clear analysis cache."""
    count = len(analysis_cache)
    analysis_cache.clear()
    return {"message": f"Cleared {count} cached analyses"}


if __name__ == "__main__":
    print("🚀 Starting MoleculeAI REST API...")
    print("📍 API running at: http://localhost:8000")
    print("📖 API docs at: http://localhost:8000/docs")
    print("⚡ Fast analysis: Only 2 API calls per molecule")
    print("")

    uvicorn.run(
        app,
        host="0.0.0.0",
        port=8000,
        log_level="info"
    )
