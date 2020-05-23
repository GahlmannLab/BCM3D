function final_matrix = rotCellWithVector(original_matrix, original_unit_vector, final_unit_vector)
% Mingxing Zhang, Gahlmann lab, Chemistry department, University of Virginia
% last edit by Mingxing 20190717
% generate the rotation vector between the original unit vector of the cell and the
% final unit vector of the cell. The rotation vector is a four-element axis-angle rotation row vector.
% The first three elements specify the rotation axis, and the last element defines the angle of rotation.
% a = [1 0 0];
% b = [0.527938425540924,0.849204242229462,0.0115396501496434]; % out put from CellModeller
% rotation_vector = vrrotvec(a, b);
rotation_vector = vrrotvec(original_unit_vector, final_unit_vector);

% convert the rotation vector to be rotation matrix. The matrix is a 3 by 3
% matrix. In order to use this matrix to rotate a cell (all points in the
% cell), the matrix should be left-multiplyed to the matrix of points in
% the cell. In addition, the format of the matrix (n rows by 3 colums) of points in the
% cell must be converted to 3 columns by n rows
rotation_matrix = vrrotvec2mat(rotation_vector);

final_matrix = rotation_matrix*original_matrix';
final_matrix = final_matrix';
end

