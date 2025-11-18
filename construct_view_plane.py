import numpy as np 
from connection_m import connection
import matplotlib.pyplot as plt 
import json 
import helper as hlp # type: ignore
import time # type: ignore

class view_plane:
    """"""
    def __init__(self, json_file):
        self.viewplane_json_file = json_file
        self.parse_viewplane_json()

    def parse_viewplane_json(self):
            """This function reads the json file and stores the cylinder specific data in a dictionary"""
            with open(self.viewplane_json_file, 'r') as json_handle:
                input = json.load(json_handle)
                # camera stuff
                self.distance_observe_to_viewplane = hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "distance[ft]"], 1.0)
                self.observation_angle = np.deg2rad(hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "angle[deg]"], 45.0))
                self.viewplane_RA = hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "aspect_ratio"], 2.0)
                self.camera_location_xyz = np.array(hlp.parse_dictionary_or_return_default(input, ["camera", "location_xyz_from_vehicle[ft]"], [0.0,0.0,0.0]))
                self.camera_orientation_phi_theta_psi = np.array(np.deg2rad(hlp.parse_dictionary_or_return_default(input, ["camera", "orientation_phi_theta_psi[deg]"], [0.0,0.0,0.0])))
                self.camera_quaternion = hlp.euler_to_quat(self.camera_orientation_phi_theta_psi)
                # ground grid
                self.ground_altitude = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "altitude[ft]"], 0.0)
                self.ground_grid_number = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "grid_number"], 0.0)
                self.ground_grid_scale = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "grid_scale[ft]"], 0.0)
                self.ground_grid_color = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "color"], 0.0)
                self.states_connection = connection(input["connections"]["receive_states"])
                self.states = [0]*14
                self.frame = 0
                self.fps = 0.0
                # aircraft
                

    def camera_set_state(self, camera_location, quat):
        """"""
        self.camera_location_xyz = camera_location
        self.camera_quaternion = quat 
        # Calculate the coordinates of the viewplane corners 
        self.calc_coordinates_of_viewplane_corners()
        # convert those to a 2D thing
        self.convert_3d_corner_coordinates_into_xy_plane()
        # calc corners in earth-fixed coordinates
        self.corners_in_earth_fixed()
        # calc PO
        self.calc_PO()
        # calc P1_P2
        self.calc_P1_P2()
        # calc_normal_vec_viewplane
        self.calc_normal_vec_viewplane()


    def calc_width_viewplane(self):
        """Calculates the width of the viewplane based on self.distance_observe_to_viewplane, self.observation_angle"""
        self.width_viewplane = 2*self.distance_observe_to_viewplane*np.tan(self.observation_angle) # Eq. 11.2.1 in flight sim book

    def calc_height_viewplane(self):
        """Calculates the height of the viewplane based on self.viewplane_RA and self.width_viewplane"""
        self.height_viewplane = self.width_viewplane/self.viewplane_RA

    def calc_coordinates_of_viewplane_corners(self):
        """Calculates the corners of the viewplane corners based on the viewplane distance away from the object, the viewplane height, and the viewplan width"""
        # create 3 4-element zero arrays for the x, y, and z coordinates of the 4 different corners of the viewplane going from top left, bottom left, bottom right, top right. (This is arbitrary, but should be kept track of)
        self.x_camera_viewplane_corners = np.array([self.distance_observe_to_viewplane, self.distance_observe_to_viewplane, self.distance_observe_to_viewplane, self.distance_observe_to_viewplane])
        self.y_camera_viewplane_corners = 0.5*np.array([-self.width_viewplane, -self.width_viewplane, self.width_viewplane, self.width_viewplane])
        self.z_camera_viewplane_corners = 0.5*np.array([-self.height_viewplane, self.height_viewplane, self.height_viewplane, -self.height_viewplane])

        # these coordinates are from Eq. 11.2.3 in the book
    
    def convert_3d_corner_coordinates_into_xy_plane(self):
        """This function converts the 3D coordinates such that they appear in an xy plane for visualization"""
        self.x_2D_corners = self.y_camera_viewplane_corners
        self.y_2D_corners = -self.z_camera_viewplane_corners
    
    def plot_viewplane_in_2D(self):
        """plots the corners in 2D"""
        lines_2D_list = []
        for i in range(self.ground_num_lines):
            first_point = self.ground_lines[i,0]
            second_point = self.ground_lines[i,1]
            if self.lambda_array[first_point] > 0 and self.lambda_array[second_point] > 0:
                self.lines2D[3*i, :] = self.projected_xy_array[first_point,:]
                self.lines2D[3*i+1,:] = self.projected_xy_array[second_point,:]
            elif self.lambda_array[first_point] < 0 and self.lambda_array[second_point] > 0:
                lca = self.ground_points[second_point]  - self.ground_points[first_point] 
                if not np.isclose(np.dot(lca, self.norm_viewplane), 0.0, atol=1e-8): # prevents division by zero
                    lamb =  np.dot(self.P0 - self.ground_points[first_point], self.norm_viewplane) / np.dot(lca, self.norm_viewplane) # this is lambda from eq. 11.3.5. The numerator has to be recomputed for the new Lca
                    pb = self.ground_points[first_point] + lamb * lca                         # 3D intersection point on plane
                    lca_pb = pb - self.camera_location_xyz   # 3D vector from camera to pb
                    lambda_pb = self.lambda_numerator/np.dot(lca_pb, self.norm_viewplane)
                    Rotation_vec = self.rotation_matrix_earth_to_body(lambda_pb, lca_pb) # 1.0 because we want the full vector from camera to pb
                    projected_point = np.array([Rotation_vec[1], -Rotation_vec[2]])
                    # set the visible end to the intersection point, the other end to the visible point's projection
                    self.lines2D[3*i, :] = projected_point
                    self.lines2D[3*i+1, :] = self.projected_xy_array[second_point, :]
            elif self.lambda_array[first_point] > 0 and self.lambda_array[second_point] < 0:
                lca = self.ground_points[second_point]  - self.ground_points[first_point] 
                if not np.isclose(np.dot(lca, self.norm_viewplane), 0.0, atol=1e-8): # prevents division by zero
                    lamb =  np.dot(self.P0 - self.ground_points[first_point], self.norm_viewplane) / np.dot(lca, self.norm_viewplane) # this is lambda from eq. 11.3.5. The numerator has to be recomputed for the new Lca
                    pb = self.ground_points[first_point] + lamb * lca # 3D intersection point on plane
                    lca_pb = pb - self.camera_location_xyz   # 3D vector from camera to pb
                    lambda_pb = self.lambda_numerator/np.dot(lca_pb, self.norm_viewplane)
                    Rotation_vec = self.rotation_matrix_earth_to_body(lambda_pb, lca_pb) # 1.0 because we want the full vector from camera to pb
                    projected_point = np.array([Rotation_vec[1], -Rotation_vec[2]])
                    # set the visible end to the visible point's projection, the other end to the intersection point
                    self.lines2D[3*i, :] = self.projected_xy_array[first_point, :]
                    self.lines2D[3*i+1, :] = projected_point
            # if one of the points is not visible, the line should extend from the visible point to the invisible point, where the point on the viewplane is pb = pc + lambda l_ca from eq. 11.3.2
        self.ax.set_data(self.lines2D[:,0], self.lines2D[:,1])

    def rotation_matrix_for_body_fixed_to_earth_fixed(self, x_cp, y_cp, z_cp):
        """"""
        vec_camera = np.array([x_cp, y_cp, z_cp])
        return hlp.quat_dependent_to_base(vec_camera, self.camera_quaternion) + self.camera_location_xyz
    
    def corners_in_earth_fixed(self): 
        """"""
        first = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[0], self.y_camera_viewplane_corners[0], self.z_camera_viewplane_corners[0])
        second = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[1], self.y_camera_viewplane_corners[1], self.z_camera_viewplane_corners[1])
        third = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[2], self.y_camera_viewplane_corners[2], self.z_camera_viewplane_corners[2])
        fourth = self.rotation_matrix_for_body_fixed_to_earth_fixed(self.x_camera_viewplane_corners[3], self.y_camera_viewplane_corners[3], self.z_camera_viewplane_corners[3])
        self.x_camera_viewplane_corners_earth_fixed = np.array([first[0], second[0], third[0], fourth[0]])
        self.y_camera_viewplane_corners_earth_fixed = np.array([first[1], second[1], third[1], fourth[1]])
        self.z_camera_viewplane_corners_earth_fixed = np.array([first[2], second[2], third[2], fourth[2]])

    def rotation_matrix_earth_to_body(self, lambda_scalar, lca_element_vector):
        """"""
        return hlp.quat_base_to_dependent(lambda_scalar*lca_element_vector,self.camera_quaternion)

    def calc_PO(self):
        """"""
        x_avg_location = np.average(self.x_camera_viewplane_corners_earth_fixed)
        y_avg_location = np.average(self.y_camera_viewplane_corners_earth_fixed)
        z_avg_location = np.average(self.z_camera_viewplane_corners_earth_fixed)
        self.P0 = np.array([x_avg_location, y_avg_location, z_avg_location])

    def calc_P1_P2(self):
        """"""
        self.P1_array = np.array([self.x_camera_viewplane_corners_earth_fixed[3],self.y_camera_viewplane_corners_earth_fixed[3],self.z_camera_viewplane_corners_earth_fixed[3]])
        self.P2_array = np.array([self.x_camera_viewplane_corners_earth_fixed[0],self.y_camera_viewplane_corners_earth_fixed[0],self.z_camera_viewplane_corners_earth_fixed[0]])

    def calc_normal_vec_viewplane(self):
        """"""
        self.norm_viewplane = np.cross(self.P1_array-self.P0, self.P2_array-self.P0)
        self.P0_minus_Pc = self.calc_length_a_minus_c(self.P0)

    def calc_length_a_minus_c(self, point_a):
        return point_a - self.camera_location_xyz
    
    def calc_l_ca_array(self, array_of_points):
        """"""
        l_ca_array_points = np.copy(array_of_points)
        for i in range(len(array_of_points)):
            l_ca_array_points[i] = self.calc_length_a_minus_c(array_of_points[i])
        return l_ca_array_points
    
    def calc_lambda_numerator(self):
        """"""
        self.lambda_numerator = np.dot(self.P0_minus_Pc, self.norm_viewplane)

    def calc_lambda_denominator(self, lca):
        """"""
        return np.dot(lca, self.norm_viewplane)
    
    def calc_lambda_array(self, array_of_lca):
        """"""
        self.calc_lambda_numerator()
        lambda_array = np.zeros(len(array_of_lca))
        for i in range(len(array_of_lca)):
            lambda_array[i] = self.lambda_numerator/(self.calc_lambda_denominator(array_of_lca[i]))
        self.lambda_array = lambda_array
        return lambda_array
    
    def calc_xy_projection_onto_viewplane(self, array_of_lca, array_lambda):
        """"""
        x_y_array = np.zeros((len(array_of_lca), 2))
        for i in range(len(array_of_lca)):
            Rotation_vec = self.rotation_matrix_earth_to_body(array_lambda[i], array_of_lca[i])
            x_y_array[i] = [Rotation_vec[1], -Rotation_vec[2]]
        self.projected_xy_array = x_y_array
        return x_y_array
    
    def calc_ground_grid(self):
        """"""
        ground_num_lines = self.ground_grid_number*2 + 1
        self.gound_n_lines = 2*ground_num_lines
        self.ground_points = np.zeros((4*ground_num_lines,3))
        self.ground_lines = np.zeros((2*ground_num_lines,2), dtype=int)
        for i in range(ground_num_lines): # x lines
            self.ground_points[2*i, :] = [-self.ground_grid_scale*self.ground_grid_number, self.ground_grid_scale*(i-self.ground_grid_number), -self.ground_altitude] # left side
            self.ground_points[2*i+1, :] = [self.ground_grid_scale*self.ground_grid_number, self.ground_grid_scale*(i-self.ground_grid_number), -self.ground_altitude] # right side
            self.ground_lines[i,0] = 2*i 
            self.ground_lines[i,1] = 2*i + 1
        for i in range(ground_num_lines): # y lines 
            self.ground_points[2*i+2*ground_num_lines,:] = [self.ground_grid_scale*(i-self.ground_grid_number), -self.ground_grid_scale*self.ground_grid_number, -self.ground_altitude] # bottom side
            self.ground_points[2*i+1+2*ground_num_lines,:] = [self.ground_grid_scale*(i-self.ground_grid_number), self.ground_grid_scale*self.ground_grid_number, -self.ground_altitude] # top side
            self.ground_lines[i+ground_num_lines,0] = 2*i + 2*ground_num_lines 
            self.ground_lines[i+ground_num_lines,1] = 2*i + 1 + 2*ground_num_lines
        self.ground_num_points = len(self.ground_points)
        self.ground_num_lines = len(self.ground_lines)
        self.lines2D = np.full((self.ground_num_lines*3,2), None, dtype =object)
    
if __name__ == "__main__":

    np.set_printoptions(formatter={'float': lambda x: f"{x:.12g}"})
    x_f_array = np.array([-10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , -10 , -5 , -5 , 0 , 0 , 5 , 5 , 10 , 10])
    y_f_array = np.array([-10 , -10 , -5 , -5 , 0 , 0 , 5 , 5 , 10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , 10])
    z_f_array = np.array([0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0])
    f_points = np.column_stack((x_f_array, y_f_array, z_f_array))
    print_stuff = False
    plot_stuff = True
    # make viewplane object (which sets self.distance_observe_to_viewplane, self.observation_angle, and self.viewplane_RA)
    viewplane_object = view_plane("construct_view_plane.json")
    viewplane_object.calc_ground_grid()
    # calculate the width of the viewplane using the distance 
    viewplane_object.calc_width_viewplane()
    # calculate height of the view plane based on the width and aspect ratio of the viewplane
    viewplane_object.calc_height_viewplane()
    # Calculate the coordinates of the viewplane corners 
    viewplane_object.camera_set_state(viewplane_object.camera_location_xyz, viewplane_object.camera_quaternion)
    # viewplane_object.calc_coordinates_of_viewplane_corners()
    l_ca_array = viewplane_object.calc_l_ca_array(viewplane_object.ground_points)
    # l_ca_array = viewplane_object.calc_l_ca_array_vectorized(f_points)
    lambda_array = viewplane_object.calc_lambda_array(l_ca_array)
    # lambda_array = viewplane_object.calc_lambda_vectorized(l_ca_array)
    ground_xy_projected_on_viewplane = viewplane_object.calc_xy_projection_onto_viewplane(l_ca_array, lambda_array)

    if print_stuff:
        print("")
        print("Distance observer to viewplane = ", viewplane_object.distance_observe_to_viewplane)
        print("Observation angle = ", np.rad2deg(viewplane_object.observation_angle))
        print("Viewplane aspect ratio = ", viewplane_object.viewplane_RA)
        print("Viewplane width = ", viewplane_object.width_viewplane)
        print("Viewplane height = ", viewplane_object.height_viewplane)
        print("x corners 3d = ", viewplane_object.x_camera_viewplane_corners)
        print("y corners 3d = ", viewplane_object.y_camera_viewplane_corners)
        print("z corners 3d = ", viewplane_object.z_camera_viewplane_corners)
        print("2D x corners = ", viewplane_object.x_2D_corners)
        print("2D y corners = ", viewplane_object.y_2D_corners)
        print("x corners earth fixed = ", viewplane_object.x_camera_viewplane_corners_earth_fixed)
        print("y corners earth fixed = ", viewplane_object.y_camera_viewplane_corners_earth_fixed)
        print("z corners earth fixed = ", viewplane_object.z_camera_viewplane_corners_earth_fixed)
        print("P0 = ", viewplane_object.P0)
        print("P1 = ", viewplane_object.P1_array)
        print("P2 = ", viewplane_object.P2_array)
        print("Normal Vec = ", viewplane_object.norm_viewplane)
        print("l_ca_array = \n", l_ca_array)
        print("lambda_array = \n", lambda_array)
        print("ground_xy_projected_on_viewplane = \n", ground_xy_projected_on_viewplane)
        print("")
    # plot everything 
    if plot_stuff:
        fig = plt.figure(figsize = (viewplane_object.viewplane_RA*5.0,5.0))
        ax = fig.add_subplot(111)
        viewplane_object.ax, = ax.plot([],[], color = viewplane_object.ground_grid_color)
        plt.subplots_adjust(top = 1.0, bottom = 0.0, left = 0.0, right = 1.0)
        plt.axis('off')
        x_corners_for_plotting = np.array([viewplane_object.x_2D_corners[0], viewplane_object.x_2D_corners[1], viewplane_object.x_2D_corners[2], viewplane_object.x_2D_corners[3], viewplane_object.x_2D_corners[0]])
        y_corners_for_plotting = np.array([viewplane_object.y_2D_corners[0], viewplane_object.y_2D_corners[1], viewplane_object.y_2D_corners[2], viewplane_object.y_2D_corners[3], viewplane_object.y_2D_corners[0]])
        ax.axes.set_xlim(x_corners_for_plotting[0], x_corners_for_plotting[2])
        ax.axes.set_ylim(y_corners_for_plotting[1], y_corners_for_plotting[0])
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axes.set_aspect('equal')
        line, = ax.plot([],[], color = viewplane_object.ground_grid_color)
        plt.show(block = False)
        while(viewplane_object.camera_location_xyz[0]<0.0):
            time_start = time.time()
            viewplane_object.plot_viewplane_in_2D()
            fig.canvas.draw()
            fig.canvas.flush_events()
            # plt.show()
            viewplane_object.camera_location_xyz[0] += 0.1
            viewplane_object.camera_set_state(viewplane_object.camera_location_xyz, viewplane_object.camera_quaternion)
            l_ca_array = viewplane_object.calc_l_ca_array(viewplane_object.ground_points)
            lambda_array = viewplane_object.calc_lambda_array(l_ca_array)
            ground_xy_projected_on_viewplane = viewplane_object.calc_xy_projection_onto_viewplane(l_ca_array, lambda_array)
            time_end = time.time()
            viewplane_object.fps = 1/(time_end-time_start)
            print("      update hz = ", viewplane_object.fps)
          