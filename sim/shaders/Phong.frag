#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  
  // (Placeholder code. You will want to replace it.)
    float ka = 0.1;
    float ks = 0.5;
    float p = 100.0;
    vec4 kd = u_color;
    vec4 Ia = vec4(1.0, 1.0, 1.0, 1.0);
    
    vec4 out_color_ambient = ka * Ia;
    
    vec3 l = normalize(u_light_pos - v_position.xyz);
    float radius_squared = distance(u_light_pos, v_position.xyz) * distance(u_light_pos, v_position.xyz);
    vec4 out_color_diffuse = kd * vec4(u_light_intensity, 0.0) / radius_squared * max(0.0, dot(normalize(v_normal.xyz), normalize(l)));
    
    vec3 v = normalize(u_cam_pos - v_position.xyz);
    vec3 h = (l + v) / length(l + v);
    vec4 out_color_specular = ks * vec4(u_light_intensity, 0.0) / radius_squared *vec4(1.0, 1.0, 1.0, 1.0) * pow(max(0.0, dot(normalize(v_normal.xyz), h)), p);
    out_color = out_color_ambient + out_color_diffuse + out_color_specular;
    out_color.a = 1;
}

