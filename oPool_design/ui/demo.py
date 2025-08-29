#!/usr/bin/env python3
"""
Demo script for testing the oPool Design Pipeline
"""

import pandas as pd
import os

def create_demo_data():
    """Create sample data for testing"""
    
    # Create demo antibody data
    demo_data = {
        'Name': ['Demo_Ab_001', 'Demo_Ab_002', 'Demo_Ab_003'],
        'VH_nuc': ['ATG...', 'ATG...', 'ATG...'],
        'VH_AA': ['QVQL...', 'QVQL...', 'QVQL...'],
        'VL_nuc': ['ATG...', 'ATG...', 'ATG...'],
        'VL_AA': ['DIQMT...', 'DIQMT...', 'DIQMT...'],
        'Heavy_V_gene': ['IGHV1-2', 'IGHV3-1', 'IGHV4-1'],
        'Heavy_J_gene': ['IGHJ1', 'IGHJ2', 'IGHJ3'],
        'Heavy_D_gene': ['IGHD1-1', 'IGHD2-1', 'IGHD3-1'],
        'Light_V_gene': ['IGLV1-1', 'IGLV2-1', 'IGLV3-1'],
        'Light_J_gene': ['IGLJ1', 'IGLJ2', 'IGLJ3'],
        'Specificity': ['HA:Unk', 'HA:Unk', 'HA:Unk']
    }
    
    df = pd.DataFrame(demo_data)
    
    # Save demo data
    os.makedirs('data', exist_ok=True)
    df.to_csv('data/demo_antibodies.csv', index=False)
    print("✅ Created demo data: data/demo_antibodies.csv")
    
    # Create demo negative control
    neg_data = {
        'Name': ['Neg_001', 'Neg_002'],
        'VH_nuc': ['ATG...', 'ATG...'],
        'VH_AA': ['QVQL...', 'QVQL...'],
        'VL_nuc': ['ATG...', 'ATG...'],
        'VL_AA': ['DIQMT...', 'DIQMT...']
    }
    
    neg_df = pd.DataFrame(neg_data)
    neg_df.to_csv('data/demo_negative.csv', index=False)
    print("✅ Created demo negative controls: data/demo_negative.csv")

if __name__ == '__main__':
    create_demo_data()
    print("\n🎉 Demo data created successfully!")
    print("You can now test the web UI with these sample files.") 