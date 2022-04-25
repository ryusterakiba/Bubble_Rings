#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
    vec4 color = texture(u_texture_2, uv);
    return color.r;
}

void main() {
  // YOUR CODE HERE
    vec4 t = v_tangent;
    vec4 n = v_normal;
    vec3 b = cross(n.xyz, t.xyz);
    
    mat3 TBN = mat3(t.xyz, b.xyz, n.xyz);
    
    float kh = u_height_scaling;
    float kn = u_normal_scaling;
    vec2 v_uv_dw = vec2(v_uv.x + 1.0 / u_texture_2_size.x, v_uv.y);
    vec2 v_uv_dh = vec2(v_uv.x, v_uv.y + 1.0 / u_texture_2_size.y);
    
    float du = (h(v_uv_dw) - h(v_uv)) * kh * kn;
    float dv = (h(v_uv_dh) - h(v_uv)) * kh * kn;
    
    vec3 localspace_normal = vec3(-du, -dv, 1.0);
    
    vec4 displaced_modelspace_normal = vec4(TBN * localspace_normal, 0.0);
    
    float ka = 0.1;
    float ks = 0.5;
    float p = 100.0;
    vec4 kd = u_color;
    vec4 Ia = vec4(1.0, 1.0, 1.0, 1.0);
    
    vec4 out_color_ambient = ka * Ia;
    
    vec3 l = normalize(u_light_pos - v_position.xyz);
    float radius_squared = distance(u_light_pos, v_position.xyz) * distance(u_light_pos, v_position.xyz);
    vec4 out_color_diffuse = kd * vec4(u_light_intensity, 0.0) / radius_squared * max(0.0, dot(normalize(displaced_modelspace_normal.xyz), normalize(l)));
    
    vec3 v = normalize(u_cam_pos - v_position.xyz);
    vec3 h = (l + v) / length(l + v);
    vec4 out_color_specular = ks * vec4(u_light_intensity, 0.0) / radius_squared *vec4(1.0, 1.0, 1.0, 1.0) * pow(max(0.0, dot(normalize(displaced_modelspace_normal.xyz), h)), p);
    //out_color = out_color_ambient + out_color_diffuse + out_color_specular;
    out_color = out_color_ambient + out_color_diffuse + out_color_specular;
    out_color.a = 1;

}

