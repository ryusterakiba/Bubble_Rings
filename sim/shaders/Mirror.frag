#version 330


uniform vec3 u_cam_pos;

uniform samplerCube u_texture_cubemap;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;

out vec4 out_color;

void main() {
    vec3 wo = u_cam_pos - v_position.xyz;
    vec3 wi = -wo + 2 * dot(normalize(wo), v_normal.xyz) * v_normal.xyz;
    
    out_color = texture(u_texture_cubemap, wi);
    out_color.a = 1;
}
