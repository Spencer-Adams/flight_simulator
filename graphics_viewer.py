import numpy as np
import sys
import pygame 
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
                # appellian stuff
                self.distance_observe_to_viewplane = hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "distance[ft]"], 1.0)
                self.observation_angle = np.deg2rad(hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "angle[deg]"], 45.0))
                self.viewplane_RA = hlp.parse_dictionary_or_return_default(input, ["camera", "view_plane", "aspect_ratio"], 2.0)
                self.camera_location_xyz = np.array(hlp.parse_dictionary_or_return_default(input, ["camera", "location_xyz_from_vehicle[ft]"], [-200.0,0.0,0.0]))
                self.original_camera_location_xyz = np.array(self.camera_location_xyz)
                self.camera_orientation_phi_theta_psi = np.array(np.deg2rad(hlp.parse_dictionary_or_return_default(input, ["camera", "orientation_phi_theta_psi[deg]"], [0.0,0.0,0.0])))
                self.camera_quaternion = hlp.euler_to_quat(self.camera_orientation_phi_theta_psi)
                # ground grid
                self.ground_altitude = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "altitude[ft]"], 0.0)
                self.ground_grid_number = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "grid_number"], 0.0)
                self.ground_grid_scale = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "grid_scale[ft]"], 0.0)
                self.ground_grid_color = hlp.parse_dictionary_or_return_default(input, ["scene", "ground", "color"], 0.0)
                # vehicle stuff 
                self.vehicle_file = hlp.parse_dictionary_or_return_default(input, ["scene", "vehicle", "vtk_file"], "F16_coarse.vtk")
                self.vehicle_location_xyz = np.array(hlp.parse_dictionary_or_return_default(input, ["scene", "vehicle", "location_xyz[ft]"], [0.0,0.0,-2000.0]))
                self.camera_location_xyz += self.vehicle_location_xyz
                print("Vehicle Location\n", self.vehicle_location_xyz)
                print("Camera Location\n", self.camera_location_xyz)
                self.ground_line = []#ax.plot([], [], color=viewplane_object.ground_grid_color)
                self.vehicle_line = []#ax.plot([], [], color='black')  # or any color
    
    def parse_vtk(self):
        filename = self.vehicle_file
        with open(filename, 'r') as f:
            lines = f.readlines()
        points_idx = None
        for i, line in enumerate(lines):
            if line.startswith("POINTS"):
                points_idx = i
                break
        if points_idx is None:
            raise RuntimeError("VTK file missing POINTS section")
        # points
        parts = lines[points_idx].split()
        n_points = int(parts[1])
        vehicle_points = np.zeros((n_points, 3))
        for i in range(n_points):
            x, y, z = map(float, lines[points_idx + 1 + i].split())
            vehicle_points[i] = [x, y, z]
        # lines 
        lines_idx = None
        for i, line in enumerate(lines):
            if line.startswith("LINES"):
                lines_idx = i
                break
        if lines_idx is None:
            raise RuntimeError("VTK file missing LINES section")
        parts = lines[lines_idx].split()
        n_lines = int(parts[1])
        vehicle_lines = np.zeros((n_lines, 2), dtype=int)
        for k in range(n_lines):
            # 2 i j
            row = lines[lines_idx + 1 + k].split()
            if int(row[0]) != 2:
                raise RuntimeError("This parser only supports 2-point LINES entries (2 i j)")
            i1 = int(row[1])
            i2 = int(row[2])
            vehicle_lines[k] = [i1, i2]
        vehicle_points += self.vehicle_location_xyz
        self.vehicle_points = vehicle_points
        # self.vehicle_points += 
        self.vehicle_lines  = vehicle_lines
        self.vehicle_num_points = n_points
        self.vehicle_num_lines  = n_lines
        # 2D projected version — same style as ground grid
        self.vehicle_lines2D = np.full((self.vehicle_num_lines * 3, 2), None, dtype=object)

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
    
    def plot_viewplane_in_2D(self, lambda_array, points, lines, num_lines, projected_xy_array, lines2D, line_type):
        """plots the corners in 2D"""
        for i in range(num_lines):
            first_point = lines[i,0]
            second_point = lines[i,1]
            if lambda_array[first_point] > 0 and lambda_array[second_point] > 0:
                lines2D[3*i, :] = projected_xy_array[first_point,:]
                lines2D[3*i+1,:] = projected_xy_array[second_point,:]
            elif lambda_array[first_point] < 0 and lambda_array[second_point] > 0:
                lca = points[second_point]  - points[first_point] 
                if not np.isclose(np.dot(lca, self.norm_viewplane), 0.0, atol=1e-8): # prevents division by zero
                    lamb =  np.dot(self.P0 - points[first_point], self.norm_viewplane) / np.dot(lca, self.norm_viewplane) # this is lambda from eq. 11.3.5. The numerator has to be recomputed for the new Lca
                    pb = points[first_point] + lamb * lca                         # 3D intersection point on plane
                    lca_pb = pb - self.camera_location_xyz   # 3D vector from camera to pb
                    lambda_pb = self.lambda_numerator/np.dot(lca_pb, self.norm_viewplane)
                    Rotation_vec = self.rotation_matrix_earth_to_body(lambda_pb, lca_pb) # 1.0 because we want the full vector from camera to pb
                    projected_point = np.array([Rotation_vec[1], -Rotation_vec[2]])
                    # set the visible end to the intersection point, the other end to the visible point's projection
                    lines2D[3*i, :] = projected_point
                    lines2D[3*i+1, :] = projected_xy_array[second_point, :]
            elif lambda_array[first_point] > 0 and lambda_array[second_point] < 0:
                lca = points[second_point]  - points[first_point] 
                if not np.isclose(np.dot(lca, self.norm_viewplane), 0.0, atol=1e-8): # prevents division by zero
                    lamb =  np.dot(self.P0 - points[first_point], self.norm_viewplane) / np.dot(lca, self.norm_viewplane) # this is lambda from eq. 11.3.5. The numerator has to be recomputed for the new Lca
                    pb = points[first_point] + lamb * lca # 3D intersection point on plane
                    lca_pb = pb - self.camera_location_xyz   # 3D vector from camera to pb
                    lambda_pb = self.lambda_numerator/np.dot(lca_pb, self.norm_viewplane)
                    Rotation_vec = self.rotation_matrix_earth_to_body(lambda_pb, lca_pb) # 1.0 because we want the full vector from camera to pb
                    projected_point = np.array([Rotation_vec[1], -Rotation_vec[2]])
                    # set the visible end to the visible point's projection, the other end to the intersection point
                    lines2D[3*i, :] = projected_xy_array[first_point, :]
                    lines2D[3*i+1, :] = projected_point
        line_type.set_data(lines2D[:,0], lines2D[:,1])
        # if overwrite_data:
        #     self.ground_line.set_data(lines2D[:,0], lines2D[:,1])
        # else:
        #     self.vehicle_line.set_data(lines2D[:,0], lines2D[:,1])
            

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
        lambda_array = np.zeros(len(array_of_lca))
        for i in range(len(array_of_lca)):
            lambda_array[i] = self.lambda_numerator/(self.calc_lambda_denominator(array_of_lca[i]))
        return lambda_array
    
    def calc_xy_projection_onto_viewplane(self, array_of_lca, array_lambda):
        """"""
        x_y_array = np.zeros((len(array_of_lca), 2))
        for i in range(len(array_of_lca)):
            Rotation_vec = self.rotation_matrix_earth_to_body(array_lambda[i], array_of_lca[i])
            x_y_array[i] = [Rotation_vec[1], -Rotation_vec[2]]
        return x_y_array

    def calc_ground_grid(self):
        scale  = self.ground_grid_scale  # spacing between lines
        N  = self.ground_grid_number # number of positive/negative steps
        Z  = -self.ground_altitude
        # calculate grid number based on ground altitude
        # if abs(Z-self.camera_location_xyz[2]) < 100.0:
        #     scale = 150.0
        # elif 100.0 <= abs(Z-self.camera_location_xyz[2]) < 2000.0:
        #     scale = 0.10101*abs(Z-self.camera_location_xyz[2]) + 150.0
        # else: 
        #     scale = 0.10101*2000.0 + 150.0
        # Aircraft current x,y
        cam_x, cam_y = self.camera_location_xyz[:2]
        # Determine which grid cell aircraft sits in
        cx = np.floor(cam_x / scale)
        cy = np.floor(cam_y / scale)
        # The new "center cell" of the grid
        center_x = cx * scale
        center_y = cy * scale
        # Range of line coordinates
        offsets = (np.arange(-N, N+1) * scale)
        # Total number of lines in each direction
        ground_num_lines = 2 * N + 1
        self.gound_n_lines = 2 * ground_num_lines
        # Preallocate (same shape as original code)
        self.ground_points = np.zeros((4 * ground_num_lines, 3))
        self.ground_lines  = np.zeros((2 * ground_num_lines, 2), dtype=int)
        # X-parallel lines (horizontal) — vary in y
        for i in range(ground_num_lines):
            y = center_y + offsets[i]
            x_left  = center_x - N * scale
            x_right = center_x + N * scale
            self.ground_points[2*i, :]   = [x_left,  y, Z]
            self.ground_points[2*i+1, :] = [x_right, y, Z]
            self.ground_lines[i, :]      = [2*i, 2*i+1]
        # Y-parallel lines (vertical) — vary in x
        base = 2 * ground_num_lines
        for i in range(ground_num_lines):
            x = center_x + offsets[i]
            y_bottom = center_y - N * scale
            y_top    = center_y + N * scale
            self.ground_points[base + 2*i,   :] = [x, y_bottom, Z]
            self.ground_points[base + 2*i+1, :] = [x, y_top,    Z]
            self.ground_lines[ground_num_lines + i, :] = [base + 2*i, base + 2*i + 1]
        # Final bookkeeping
        self.ground_num_points = len(self.ground_points)
        self.ground_num_lines  = len(self.ground_lines)
        self.lines2D = np.full((self.ground_num_lines * 3, 2), None, dtype=object)
    
if __name__ == "__main__":

    np.set_printoptions(formatter={'float': lambda x: f"{x:.12g}"})
    # x_f_array = np.array([-10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , -10 , -5 , -5 , 0 , 0 , 5 , 5 , 10 , 10])
    # y_f_array = np.array([-10 , -10 , -5 , -5 , 0 , 0 , 5 , 5 , 10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , 10 , -10 , 10])
    # z_f_array = np.array([0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0])
    # f_points = np.column_stack((x_f_array, y_f_array, z_f_array))
    print_stuff = False
    plot_stuff = True
    # make viewplane object (which sets self.distance_observe_to_viewplane, self.observation_angle, and self.viewplane_RA)
    viewplane_object = view_plane("graphics.json")

    with open(viewplane_object.viewplane_json_file, 'r') as json_handle:
        file_loc = json.load(json_handle)
    states_connection = connection(file_loc["connections"]["receive_states"])
    controls_connection = connection(file_loc["connections"]["send_states"])
    pygame.init()
    pygame.joystick.init()
    joy = pygame.joystick.Joystick(0)
    joy.init()
    ## axis 0: Rudder
    ## axis 1: Throttle
    ## axis 2: Aileron 
    ## axis 3: elevator
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            pygame.quit()
            sys.exit() 
    states = [0]*14
    frame = 0
    fps = 0.0
    # parse vehicle vtk
    viewplane_object.parse_vtk()
    # calculate the width of the viewplane using the distance 
    viewplane_object.calc_width_viewplane()
    # calculate height of the view plane based on the width and aspect ratio of the viewplane
    viewplane_object.calc_height_viewplane()
    viewplane_object.camera_set_state(viewplane_object.camera_location_xyz, viewplane_object.camera_quaternion)
    viewplane_object.calc_ground_grid()
    l_ca_array = viewplane_object.calc_l_ca_array(viewplane_object.ground_points)
    l_ca_vehicle_array = viewplane_object.calc_l_ca_array(viewplane_object.vehicle_points)
    viewplane_object.calc_lambda_numerator()
    lambda_array = viewplane_object.calc_lambda_array(l_ca_array)
    lambda_vehicle_array = viewplane_object.calc_lambda_array(l_ca_vehicle_array)
    ground_xy_projected_on_viewplane = viewplane_object.calc_xy_projection_onto_viewplane(l_ca_array, lambda_array)
    vehicle_xy_projected_on_viewplane = viewplane_object.calc_xy_projection_onto_viewplane(l_ca_vehicle_array, lambda_vehicle_array)

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
    if plot_stuff:
        fig = plt.figure(figsize=(viewplane_object.viewplane_RA*5.0, 5.0))
        ax = fig.add_subplot(111)
        # Store the Axes object
        viewplane_object.ax = ax
        plt.subplots_adjust(top=1, bottom=0, left=0, right=1)
        plt.axis('off')
        # Create separate line handles
        viewplane_object.ground_line, = ax.plot([], [], color=viewplane_object.ground_grid_color)
        viewplane_object.vehicle_line, = ax.plot([], [], color='black')  # or any color
        ax.set_xlim(viewplane_object.x_2D_corners[0], viewplane_object.x_2D_corners[2])
        ax.set_ylim(viewplane_object.y_2D_corners[1], viewplane_object.y_2D_corners[0])
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_aspect('equal')
        # Draw ground
        plt.show(block = False)
        viewplane_object.plot_viewplane_in_2D(lambda_array,viewplane_object.ground_points,viewplane_object.ground_lines,viewplane_object.ground_num_lines,ground_xy_projected_on_viewplane,viewplane_object.lines2D,line_type = viewplane_object.ground_line)
        # Draw vehicle
        viewplane_object.plot_viewplane_in_2D(lambda_vehicle_array,viewplane_object.vehicle_points,viewplane_object.vehicle_lines,viewplane_object.vehicle_num_lines,vehicle_xy_projected_on_viewplane,viewplane_object.vehicle_lines2D,line_type = viewplane_object.vehicle_line)
        # plt.show()
        while(frame<1000.0):
            time_start = time.time()
            viewplane_object.plot_viewplane_in_2D(lambda_array, viewplane_object.ground_points, viewplane_object.ground_lines, viewplane_object.ground_num_lines, ground_xy_projected_on_viewplane, viewplane_object.lines2D, viewplane_object.ground_line)
            # viewplane_object.plot_viewplane_in_2D(lambda_vehicle_array, viewplane_object.vehicle_points, viewplane_object.vehicle_lines, viewplane_object.vehicle_num_lines, vehicle_xy_projected_on_viewplane, viewplane_object.vehicle_lines2D)
            fig.canvas.draw()
            fig.canvas.flush_events()
            Controls = np.array([
                np.round(-joy.get_axis(2)*0.375246,1),  # da -> val * 21.5 * Pi/180
                np.round(-joy.get_axis(3)*0.436332,1),  # de -> val * 25.0 * Pi/180
                np.round(joy.get_axis(0)*0.523599,1),  # dr -> val * 30.0 * Pi/180
                np.round(-0.5*joy.get_axis(1)+0.5,1),  # tau -> val * -0.5 + 0.5 so that joystick goes from 0 to 1 when deflected from bottom to top
            ])
            Controls_sent = controls_connection.send(Controls)
            states = states_connection.recv()
            # print(states)
            # viewplane_object.vehicle_location_xyz[:] = states[7:10]
            viewplane_object.camera_location_xyz[:] = states[7:10] #+ viewplane_object.original_camera_location_xyz
            viewplane_object.camera_quaternion[:] = states[10:14]
            # viewplane_object.camera_location_xyz[0] += 0.1
            viewplane_object.camera_set_state(viewplane_object.camera_location_xyz, viewplane_object.camera_quaternion)
            viewplane_object.calc_ground_grid()
            l_ca_array = viewplane_object.calc_l_ca_array(viewplane_object.ground_points)
            # l_ca_vehicle_array = viewplane_object.calc_l_ca_array(viewplane_object.vehicle_points)
            viewplane_object.calc_lambda_numerator()
            lambda_array = viewplane_object.calc_lambda_array(l_ca_array)
            # lambda_vehicle_array = viewplane_object.calc_lambda_array(l_ca_vehicle_array)
            ground_xy_projected_on_viewplane = viewplane_object.calc_xy_projection_onto_viewplane(l_ca_array, lambda_array)
            time_end = time.time()
            fps = 1/(time_end-time_start)
            print("      update hz = ", fps)
            # print("        aileron = ", Controls[0]*180/np.pi, " elevator = ", Controls[1]*180/np.pi, " rudder = ", Controls[2]*180/np.pi, " throttle = ", Controls[3], "update hz = ", fps)
