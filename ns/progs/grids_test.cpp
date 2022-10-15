#include "common.hpp"
#include "prog_common.hpp"
#include "grid/grid.hpp"
#include "grid/grid_builder.hpp"

bool equal_cell_face(const std::vector<Grid::CellFaceEntry>& cf,
		const std::vector<int>& faces,
		const std::vector<int>& directs){
	
	if (cf.size() != faces.size()) return false;
	if (cf.size() != directs.size()) return false;

	for (auto it: cf){
		int iface1 = it.face_index;
		int idirect1 = it.normal_direction;
		auto fnd = std::find(faces.begin(), faces.end(), iface1);
		if (fnd == faces.end()) return false;
		int index = fnd - faces.begin();
		int iface2 = faces[index];
		int idirect2 = directs[index];

		if (iface1 != iface2) return false;
		if (idirect1 != idirect2) return false;
	}

	return true;
}

bool equal_face_cell(const Grid::FaceCellEntry& fc, int left, int right){
	return fc.left_cell == left && fc.right_cell == right;
}


void check_gmshvtk_2d(){
	// grid
	std::string grid_filename = from_input_path("rect1.vtk");
	std::shared_ptr<Grid> grid = GridBuilder::build_from_gmshvtk(grid_filename);

	// check sizes
	CHECK(grid->dim == 2);
	CHECK(grid->n_points() == 655);
	CHECK(grid->n_cells() == 1156);
	CHECK(grid->n_faces() == 1810);

	// check connectivities
	// cell->point
	CHECK(equal_int_vec(grid->tab_cell_point(0), {401, 624, 554}));
	CHECK(equal_int_vec(grid->tab_cell_point(300), {443, 444, 171}));
	// point->cell
	CHECK(equal_int_vec_anyorder(grid->tab_point_cell(0), {891, 1136}));
	CHECK(equal_int_vec_anyorder(grid->tab_point_cell(300), {181, 860, 753, 710, 811, 806}));
	// face->point
	CHECK(equal_int_vec_anyorder(grid->tab_face_point(0), {0, 4}));
	CHECK(equal_int_vec_anyorder(grid->tab_face_point(300), {99, 563}));
	// point->face
	CHECK(equal_int_vec_anyorder(grid->tab_point_face(0), {0, 1, 2}));
	CHECK(equal_int_vec_anyorder(grid->tab_point_face(300), {728, 732, 1062, 1059, 1061, 1060}));
	// cell->face
	CHECK(equal_cell_face(grid->tab_cell_face(0), {1454, 1755, 1453}, {1, -1, -1})); 
	CHECK(equal_cell_face(grid->tab_cell_face(300), {537, 536, 1564}, {-1, 1, 1})); 
	// face->cell
	CHECK(equal_face_cell(grid->tab_face_cell(0), 891, -1))
	CHECK(equal_face_cell(grid->tab_face_cell(1325), 1009, 239))

	// check boundary
	CHECK(equal_int_vec_anyorder(grid->btypes(), {1, 2, 3, 4}));
	auto b1 = grid->boundary(1);
	auto tf = b1.tab_faces();
	CHECK(b1.n_faces() == 8);
	CHECK(equal_int_vec_anyorder(tf, {1, 451, 448, 445, 442, 439, 436, 10}));
	CHECK(equal_int_vec_anyorder(b1.tab_points(), {3, 145, 146, 147, 148, 149, 150, 151, 0}));
	int f1 = std::find(tf.begin(), tf.end(), 1) - tf.begin();
	CHECK(equal_int_vec(b1.tab_face_point_positive(f1), {151, 0}));
	int f442 = std::find(tf.begin(), tf.end(), 442) - tf.begin();
	CHECK(equal_int_vec(b1.tab_face_point_positive(f442), {147, 148}));

	// save_vtk
	grid->save_vtk(from_output_path("g.vtk"));
	grid->save_vtk_faces(from_output_path("g_faces.vtk"));
	
	// read from saved
	std::shared_ptr<Grid> grid2 = GridBuilder::build_from_gmshvtk(from_output_path("g.vtk"));
	CHECK(grid2->dim == 2);
	CHECK(grid2->n_points() == 655);
	CHECK(grid2->n_cells() == 1156);
	CHECK(grid2->n_faces() == 1810);
	CHECK(grid2->btypes().empty());
}

int main(){
	try{
		check_gmshvtk_2d();

		std::cout << "DONE" << std::endl;
	} catch (std::exception& e){
		std::cout << "ERROR: " << "  " << e.what() << std::endl;
	}
}
