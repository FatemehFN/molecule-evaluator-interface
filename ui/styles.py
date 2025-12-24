"""Custom CSS styles for the application."""
import streamlit as st


def apply_custom_styles():
    """Apply custom CSS styles to the Streamlit app."""
    st.markdown("""
    <style>
        /* Main container */
        .main {
            padding: 2rem;
        }

        /* Header styling */
        .header-container {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 2rem;
            border-radius: 10px;
            color: white;
            margin-bottom: 2rem;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }

        .header-title {
            font-size: 2.5rem;
            font-weight: 700;
            margin: 0;
            text-align: center;
        }

        .header-subtitle {
            font-size: 1.2rem;
            text-align: center;
            margin-top: 0.5rem;
            opacity: 0.9;
        }

        /* Molecule card styling */
        .molecule-card {
            background: white;
            border: 1px solid #e0e0e0;
            border-radius: 12px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.05);
            transition: transform 0.2s, box-shadow 0.2s;
        }

        .molecule-card:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.1);
        }

        .molecule-name {
            font-size: 1.5rem;
            font-weight: 600;
            color: #333;
            margin-bottom: 1rem;
        }

        .property-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 1rem;
            margin: 1rem 0;
        }

        .property-item {
            background: #f8f9fa;
            padding: 0.75rem;
            border-radius: 6px;
            border-left: 3px solid #667eea;
        }

        .property-label {
            font-size: 0.85rem;
            color: #666;
            font-weight: 500;
            margin-bottom: 0.25rem;
        }

        .property-value {
            font-size: 1rem;
            color: #333;
            font-weight: 600;
        }

        /* Chat interface styling */
        .chat-container {
            background: #f8f9fa;
            border-radius: 10px;
            padding: 1rem;
            margin: 1rem 0;
            max-height: 600px;
            overflow-y: auto;
        }

        .chat-message {
            padding: 1rem;
            margin: 0.5rem 0;
            border-radius: 8px;
            max-width: 80%;
        }

        .user-message {
            background: #667eea;
            color: white;
            margin-left: auto;
            text-align: right;
        }

        .assistant-message {
            background: white;
            border: 1px solid #e0e0e0;
            margin-right: auto;
        }

        /* Section styling */
        .section-header {
            font-size: 1.3rem;
            font-weight: 600;
            color: #333;
            margin: 1.5rem 0 1rem 0;
            padding-bottom: 0.5rem;
            border-bottom: 2px solid #667eea;
        }

        /* Info boxes */
        .info-box {
            background: #e3f2fd;
            border-left: 4px solid #2196f3;
            padding: 1rem;
            border-radius: 4px;
            margin: 1rem 0;
        }

        .success-box {
            background: #e8f5e9;
            border-left: 4px solid #4caf50;
            padding: 1rem;
            border-radius: 4px;
            margin: 1rem 0;
        }

        .warning-box {
            background: #fff3e0;
            border-left: 4px solid #ff9800;
            padding: 1rem;
            border-radius: 4px;
            margin: 1rem 0;
        }

        .error-box {
            background: #ffebee;
            border-left: 4px solid #f44336;
            padding: 1rem;
            border-radius: 4px;
            margin: 1rem 0;
        }

        /* Button styling */
        .stButton > button {
            border-radius: 8px;
            font-weight: 500;
            transition: all 0.3s;
        }

        .stButton > button:hover {
            transform: translateY(-1px);
            box-shadow: 0 4px 8px rgba(0, 0, 0, 0.15);
        }

        /* Metrics styling */
        [data-testid="stMetricValue"] {
            font-size: 1.8rem;
            font-weight: 700;
        }

        /* Tabs styling */
        .stTabs [data-baseweb="tab-list"] {
            gap: 8px;
        }

        .stTabs [data-baseweb="tab"] {
            border-radius: 8px 8px 0 0;
            padding: 0.5rem 1.5rem;
            font-weight: 500;
        }

        /* Sidebar styling */
        [data-testid="stSidebar"] {
            background: linear-gradient(180deg, #f8f9fa 0%, #ffffff 100%);
        }

        /* Code blocks */
        code {
            background: #f4f4f4;
            padding: 0.2rem 0.4rem;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }

        /* Scrollbar styling */
        ::-webkit-scrollbar {
            width: 8px;
            height: 8px;
        }

        ::-webkit-scrollbar-track {
            background: #f1f1f1;
            border-radius: 10px;
        }

        ::-webkit-scrollbar-thumb {
            background: #888;
            border-radius: 10px;
        }

        ::-webkit-scrollbar-thumb:hover {
            background: #555;
        }
    </style>
    """, unsafe_allow_html=True)
