import text_manipulation

text_manipulation.check_lines_equal_error("mesh_spheres.geo", "ref_mesh_spheres.geo", start_line=582, end_line=2818)
text_manipulation.check_lines_equal_error("mesh_spheres_OCC.geo", "ref_mesh_spheres_OCC.geo", start_line=582, end_line=2818)
