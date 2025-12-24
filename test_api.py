#!/usr/bin/env python
"""Test suite for API connectivity and agent functionality."""

import sys
import time
import json
from config import Config


def test_1_config_load():
    """Test 1: Can we load config and API key?"""
    print("\n" + "=" * 60)
    print("TEST 1: Config Loading")
    print("=" * 60)
    
    try:
        Config.validate()
        key = Config.GEMINI_API_KEY
        if key and len(key) > 10:
            print(f"✅ Config loaded successfully")
            print(f"   API Key (first 20 chars): {key[:20]}...")
            return True
        else:
            print(f"❌ API key not set or invalid")
            return False
    except Exception as e:
        print(f"❌ Config validation failed: {e}")
        return False


def test_2_genai_import():
    """Test 2: Can we import and configure genai?"""
    print("\n" + "=" * 60)
    print("TEST 2: Genai Import & Configuration")
    print("=" * 60)
    
    try:
        import google.generativeai as genai
        print(f"✅ google.generativeai imported")
        
        Config.validate()
        genai.configure(api_key=Config.GEMINI_API_KEY)
        print(f"✅ genai.configure() successful")
        
        return True
    except Exception as e:
        print(f"❌ Failed: {e}")
        return False


def test_3_gemini_api_call():
    """Test 3: Can we make a simple API call?"""
    print("\n" + "=" * 60)
    print("TEST 3: Simple Gemini API Call (10 seconds max)")
    print("=" * 60)
    
    try:
        import google.generativeai as genai
        
        Config.validate()
        genai.configure(api_key=Config.GEMINI_API_KEY)
        
        model = genai.GenerativeModel(
            model_name=Config.GEMINI_MODEL,
            generation_config={
                "temperature": Config.TEMPERATURE,
                "max_output_tokens": 100,
            }
        )
        
        print(f"⏳ Calling Gemini API with simple prompt...")
        start = time.time()
        response = model.generate_content("Reply with just 'OK'")
        elapsed = time.time() - start
        
        if response and response.text:
            print(f"✅ API responded in {elapsed:.1f} seconds")
            print(f"   Response: {response.text[:100]}")
            return True
        else:
            print(f"❌ API returned empty response")
            return False
            
    except Exception as e:
        print(f"❌ API call failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_4_json_extraction():
    """Test 4: Can the agent extract JSON from responses?"""
    print("\n" + "=" * 60)
    print("TEST 4: JSON Extraction")
    print("=" * 60)
    
    try:
        import google.generativeai as genai
        
        Config.validate()
        genai.configure(api_key=Config.GEMINI_API_KEY)
        
        from agents.base_agent import BaseAgent
        
        class TestAgent(BaseAgent):
            def process(self):
                pass
        
        agent = TestAgent("test", "test")
        
        # Test JSON extraction
        test_response = '{"key": "value", "number": 123}'
        extracted = agent._extract_json(test_response)
        
        if extracted and extracted.get("key") == "value":
            print(f"✅ JSON extraction works")
            return True
        else:
            print(f"❌ JSON extraction failed")
            return False
            
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_5_rdkit_properties():
    """Test 5: Can RDKit calculate molecular properties?"""
    print("\n" + "=" * 60)
    print("TEST 5: RDKit Molecular Properties")
    print("=" * 60)
    
    try:
        from agents.unified_agent_fast import UnifiedAgentFast
        
        agent = UnifiedAgentFast()
        
        # Test with a known SMILES
        smiles = "CC(=O)Oc1ccccc1C(=O)O"  # Aspirin
        props = agent.get_molecular_properties_local(smiles)
        
        if props and props.get("Molecular Weight"):
            print(f"✅ RDKit working")
            print(f"   Aspirin MW: {props.get('Molecular Weight')} g/mol")
            print(f"   Formula: {props.get('Molecular Formula')}")
            return True
        else:
            print(f"❌ RDKit calculation failed")
            return False
            
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_6_cached_molecule():
    """Test 6: Does cached data work?"""
    print("\n" + "=" * 60)
    print("TEST 6: Cached Molecule (Aspirin)")
    print("=" * 60)
    
    try:
        import google.generativeai as genai
        
        Config.validate()
        genai.configure(api_key=Config.GEMINI_API_KEY)
        
        from agents.unified_agent_fast import UnifiedAgentFast
        
        agent = UnifiedAgentFast()
        
        print(f"⏳ Processing Aspirin (should be instant from cache)...")
        start = time.time()
        result = agent.process("Aspirin")
        elapsed = time.time() - start
        
        if result.get("status") == "success":
            print(f"✅ Cached molecule works ({elapsed:.2f} seconds)")
            print(f"   SMILES: {result.get('smiles')}")
            print(f"   MW: {result.get('properties', {}).get('Molecular Weight')} g/mol")
            print(f"   Cost: {result.get('market', {}).get('total_price_toman')} Toman")
            return True
        else:
            print(f"❌ Failed: {result.get('error')}")
            return False
            
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_7_unknown_molecule_api():
    """Test 7: App behavior for unknown molecules."""
    print("\n" + "=" * 60)
    print("TEST 7: Unknown Molecule Handling")
    print("=" * 60)
    print("Note: Since Gemini API is unresponsive, app returns cached-only message")
    
    try:
        import google.generativeai as genai
        
        Config.validate()
        genai.configure(api_key=Config.GEMINI_API_KEY)
        
        from agents.unified_agent_fast import UnifiedAgentFast
        
        agent = UnifiedAgentFast()
        
        print(f"⏳ Testing with unknown molecule 'Caffeine'...")
        result = agent.process("Caffeine")
        
        if result.get("status") == "unavailable":
            print(f"✅ Gracefully handled unknown molecule")
            print(f"   Message: {result.get('error')}")
            return True
        else:
            print(f"⚠️  Unexpected response: {result.get('status')}")
            return False
            
    except Exception as e:
        print(f"❌ Failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all tests and report."""
    print("\n" + "🧪 MOLECULE EVALUATOR API TEST SUITE 🧪".center(60))
    
    tests = [
        ("Config Loading", test_1_config_load),
        ("Genai Import", test_2_genai_import),
        ("Gemini API Call", test_3_gemini_api_call),
        ("JSON Extraction", test_4_json_extraction),
        ("RDKit Properties", test_5_rdkit_properties),
        ("Cached Molecule", test_6_cached_molecule),
        ("Unknown Molecule API", test_7_unknown_molecule_api),
    ]
    
    results = {}
    for name, test_func in tests:
        try:
            results[name] = test_func()
        except Exception as e:
            print(f"\n❌ Test crashed: {e}")
            results[name] = False
    
    # Summary
    print("\n" + "=" * 60)
    print("TEST SUMMARY")
    print("=" * 60)
    
    passed = sum(1 for v in results.values() if v)
    total = len(results)
    
    for name, result in results.items():
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! App is ready.")
        return 0
    elif passed >= 6:
        print(f"\n⚠️  {total - passed} test(s) failed. App may have limited functionality.")
        return 1
    else:
        print(f"\n❌ Critical tests failed. Do not use app yet.")
        return 2


if __name__ == "__main__":
    sys.exit(main())
