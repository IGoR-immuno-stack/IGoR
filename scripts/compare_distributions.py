#!/usr/bin/env python3
"""
Compare realization distributions between Original and MaxSpeed generation.
"""

import sys
import os
from collections import Counter
import statistics

def parse_realizations(filepath):
    """Parse realizations file to extract key values."""
    data = {
        'v_del': [],
        'd_del_5': [],
        'd_del_3': [],
        'j_del': [],
        'vd_ins': [],
        'dj_ins': [],
    }
    
    with open(filepath, 'r') as f:
        header = next(f).strip().split(';')
        # Find column indices
        col_map = {name: i for i, name in enumerate(header)}
        
        for line in f:
            parts = line.strip().split(';')
            if len(parts) < 8:
                continue
            
            # Extract values - they're in format "(N)" or "(N,N,N,...)"
            def extract_first_int(s):
                s = s.strip('()')
                if not s:
                    return None
                vals = s.split(',')
                return int(vals[0]) if vals[0] else None
            
            def extract_length(s):
                """Get length of comma-separated list in parentheses."""
                s = s.strip('()')
                if not s:
                    return 0
                return len(s.split(','))
            
            # Column order: seq_index, V_gene, J_gene, D_gene, V_del, D_del5, D_del3, J_del, VD_ins, VD_dinuc, DJ_ins, DJ_dinuc, Errors
            # Indices:       0         1       2       3       4      5       6       7      8       9         10      11        12
            
            v_del = extract_first_int(parts[4])
            d_del_5 = extract_first_int(parts[5])
            d_del_3 = extract_first_int(parts[6])
            j_del = extract_first_int(parts[7])
            vd_ins = extract_first_int(parts[8])
            dj_ins = extract_first_int(parts[10])
            
            if v_del is not None:
                data['v_del'].append(v_del)
            if d_del_5 is not None:
                data['d_del_5'].append(d_del_5)
            if d_del_3 is not None:
                data['d_del_3'].append(d_del_3)
            if j_del is not None:
                data['j_del'].append(j_del)
            if vd_ins is not None:
                data['vd_ins'].append(vd_ins)
            if dj_ins is not None:
                data['dj_ins'].append(dj_ins)
    
    return data


def print_comparison(name, orig_vals, max_vals):
    """Print comparison of two distributions."""
    if not orig_vals or not max_vals:
        print(f"  {name}: No data")
        return
    
    orig_mean = statistics.mean(orig_vals)
    max_mean = statistics.mean(max_vals)
    diff = max_mean - orig_mean
    pct = 100 * diff / orig_mean if orig_mean > 0 else 0
    
    print(f"  {name:12s} Orig mean={orig_mean:6.2f}, Max mean={max_mean:6.2f}, diff={diff:+.2f} ({pct:+.1f}%)")
    
    # Show top values
    orig_top = Counter(orig_vals).most_common(5)
    max_top = Counter(max_vals).most_common(5)
    print(f"              Orig top: {orig_top}")
    print(f"              Max  top: {max_top}")


def main():
    base_dir = "/tmp/igor_benchmark"
    
    orig_file = f"{base_dir}/original_1000_generated/generated_realizations_werr.csv"
    maxspeed_file = f"{base_dir}/maxspeed_1t_1000_generated/generated_realizations_werr.csv"
    
    if not os.path.exists(orig_file):
        print(f"Original file not found: {orig_file}")
        return 1
    if not os.path.exists(maxspeed_file):
        print(f"MaxSpeed file not found: {maxspeed_file}")
        return 1
    
    orig_data = parse_realizations(orig_file)
    max_data = parse_realizations(maxspeed_file)
    
    print("Distribution Comparison: Original vs MaxSpeed (1000 sequences)")
    print("=" * 70)
    
    for key in ['v_del', 'd_del_5', 'd_del_3', 'j_del', 'vd_ins', 'dj_ins']:
        print_comparison(key, orig_data[key], max_data[key])
        print()
    
    # Overall sequence length contribution from insertions and deletions
    # Length = len(V) - v_del + vd_ins + len(D) - d_del_5 - d_del_3 + dj_ins + len(J) - j_del
    # The difference in total is: -(v_del_diff) + vd_ins_diff - d_del_5_diff - d_del_3_diff + dj_ins_diff - j_del_diff
    
    print("=" * 70)
    print("Net effect on sequence length:")
    
    total_orig = sum(orig_data['vd_ins']) + sum(orig_data['dj_ins']) - sum(orig_data['v_del']) - sum(orig_data['d_del_5']) - sum(orig_data['d_del_3']) - sum(orig_data['j_del'])
    total_max = sum(max_data['vd_ins']) + sum(max_data['dj_ins']) - sum(max_data['v_del']) - sum(max_data['d_del_5']) - sum(max_data['d_del_3']) - sum(max_data['j_del'])
    
    n = len(orig_data['v_del'])
    print(f"  Original total net (ins - del): {total_orig:+d} ({total_orig/n:+.1f} per seq)")
    print(f"  MaxSpeed total net (ins - del): {total_max:+d} ({total_max/n:+.1f} per seq)")
    print(f"  Difference: {total_max - total_orig:+d} ({(total_max - total_orig)/n:+.1f} per seq)")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
